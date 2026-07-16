#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import math
import random
import re
from collections import defaultdict
from pathlib import Path


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def clean_run_column(column: str) -> str:
    column = re.sub(r"_(Abundance|Coverage)(-RPKs)?$", "", column)
    column = re.sub(r"\.(pathabundance|pathcoverage)$", "", column)
    return column


def mean(xs: list[float]) -> float:
    return sum(xs) / len(xs) if xs else float("nan")


def median(xs: list[float]) -> float:
    if not xs:
        return float("nan")
    vals = sorted(xs)
    mid = len(vals) // 2
    if len(vals) % 2:
        return vals[mid]
    return (vals[mid - 1] + vals[mid]) / 2


def cliffs_delta(xs: list[float], ys: list[float]) -> float:
    if not xs or not ys:
        return float("nan")
    gt = 0
    lt = 0
    for x in xs:
        for y in ys:
            if y > x:
                gt += 1
            elif y < x:
                lt += 1
    return (gt - lt) / (len(xs) * len(ys))


def exact_permutation_p(values: list[float], groups: list[str]) -> float:
    s_count = sum(1 for g in groups if g == "S")
    m_count = len(groups) - s_count
    if s_count == 0 or m_count == 0:
        return float("nan")
    s_idx = [i for i, g in enumerate(groups) if g == "S"]
    obs_s = sum(values[i] for i in s_idx) / s_count
    obs_m = (sum(values) - sum(values[i] for i in s_idx)) / m_count
    obs = obs_s - obs_m
    extreme = 0
    total = 0
    all_sum = sum(values)
    for subset in itertools.combinations(range(len(values)), s_count):
        s_sum = sum(values[i] for i in subset)
        diff = (s_sum / s_count) - ((all_sum - s_sum) / m_count)
        if abs(diff) >= abs(obs) - 1e-15:
            extreme += 1
        total += 1
    return extreme / total


def bootstrap_ci(xs: list[float], ys: list[float], seed: int, n_boot: int) -> tuple[float, float]:
    if not xs or not ys or n_boot <= 0:
        return (float("nan"), float("nan"))
    rng = random.Random(seed)
    diffs: list[float] = []
    for _ in range(n_boot):
        xb = [xs[rng.randrange(len(xs))] for _ in xs]
        yb = [ys[rng.randrange(len(ys))] for _ in ys]
        diffs.append(mean(yb) - mean(xb))
    diffs.sort()
    lo = diffs[int(0.025 * (len(diffs) - 1))]
    hi = diffs[int(0.975 * (len(diffs) - 1))]
    return lo, hi


def read_modules(path: Path, membership: str) -> dict[str, set[str]]:
    rows = read_tsv(path)
    modules: dict[str, set[str]] = defaultdict(set)
    accepted_status = {"finalised", "draft_seed_not_final", "locked", "draft_seed_not_final_locked"}
    for row in rows:
        if row.get("lock_status") not in accepted_status:
            continue
        if row.get("membership") != membership:
            continue
        pathway = row["pathway_id"].strip()
        module = row["module"].strip()
        if not pathway or not module:
            continue
        modules[module].add(pathway)
    fermentation: set[str] = set()
    for module, features in modules.items():
        if module.startswith("scfa_"):
            fermentation.update(features)
    if fermentation:
        modules["overall_fermentation"] = fermentation
    return dict(sorted(modules.items()))


def read_pathway_abundance(path: Path) -> tuple[list[str], dict[str, dict[str, float]]]:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        feature_col = header[0]
        if feature_col not in {"# Pathway", "feature"}:
            raise SystemExit(f"Unexpected feature column: {feature_col}")
        run_cols = [clean_run_column(c) for c in header[1:]]
        values: dict[str, dict[str, float]] = {}
        totals = {run: 0.0 for run in run_cols}
        for parts in reader:
            if not parts:
                continue
            feature = parts[0]
            if "|" in feature or feature in {"UNMAPPED", "UNINTEGRATED"}:
                continue
            row: dict[str, float] = {}
            for run, raw in zip(run_cols, parts[1:]):
                try:
                    val = float(raw)
                except ValueError:
                    val = 0.0
                row[run] = val
                totals[run] += val
            values[feature] = row
    rel: dict[str, dict[str, float]] = {}
    for feature, row in values.items():
        rel[feature] = {}
        for run, val in row.items():
            total = totals.get(run, 0.0)
            rel[feature][run] = val / total if total > 0 else 0.0
    return run_cols, rel


def resolve_module_features(module_features: set[str], available_features: set[str]) -> list[str]:
    by_id: dict[str, str] = {}
    for feature in sorted(available_features):
        pathway_id = feature.split(": ", 1)[0]
        by_id[pathway_id] = feature
    resolved: list[str] = []
    for feature in sorted(module_features):
        if feature in available_features:
            resolved.append(feature)
        elif feature in by_id:
            resolved.append(by_id[feature])
    return resolved


def build_scores(
    runs: list[str],
    rel: dict[str, dict[str, float]],
    modules: dict[str, set[str]],
    groups: dict[str, str],
    pseudocount: float,
) -> tuple[list[dict[str, str]], dict[str, dict[str, float]], dict[str, int]]:
    score_rows: list[dict[str, str]] = []
    values_by_module: dict[str, dict[str, float]] = {}
    present_counts: dict[str, int] = {}
    available_features = set(rel)
    for module, features in modules.items():
        present = resolve_module_features(features, available_features)
        present_counts[module] = len(present)
        values_by_module[module] = {}
        for run in runs:
            if groups.get(run) not in {"M", "S"}:
                continue
            raw_score = sum(rel[feature].get(run, 0.0) for feature in present)
            background = max(1.0 - raw_score, 0.0)
            log_ratio = math.log(raw_score + pseudocount) - math.log(background + pseudocount)
            values_by_module[module][run] = log_ratio
            score_rows.append(
                {
                    "run_accession": run,
                    "group": groups[run],
                    "module": module,
                    "pathways_defined": str(len(features)),
                    "pathways_present": str(len(present)),
                    "relative_score": f"{raw_score:.12g}",
                    "log_ratio_score": f"{log_ratio:.12g}",
                }
            )
    return score_rows, values_by_module, present_counts


def score_one_module(
    module_features: set[str],
    runs: list[str],
    rel: dict[str, dict[str, float]],
    pseudocount: float,
) -> dict[str, float]:
    present = resolve_module_features(module_features, set(rel))
    out: dict[str, float] = {}
    for run in runs:
        raw_score = sum(rel[feature].get(run, 0.0) for feature in present)
        background = max(1.0 - raw_score, 0.0)
        out[run] = math.log(raw_score + pseudocount) - math.log(background + pseudocount)
    return out


def module_stats(
    values_by_module: dict[str, dict[str, float]],
    groups: dict[str, str],
    module_sizes: dict[str, int],
    modules: dict[str, set[str]],
    seed: int,
    n_boot: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    stat_rows: list[dict[str, str]] = []
    loo_rows: list[dict[str, str]] = []
    for module, run_values in values_by_module.items():
        runs = sorted(run_values)
        vals = [run_values[r] for r in runs]
        grps = [groups[r] for r in runs]
        m_vals = [run_values[r] for r in runs if groups[r] == "M"]
        s_vals = [run_values[r] for r in runs if groups[r] == "S"]
        delta = mean(s_vals) - mean(m_vals)
        p_exact = exact_permutation_p(vals, grps)
        ci_lo, ci_hi = bootstrap_ci(m_vals, s_vals, seed + len(module), n_boot)
        full_sign = 0 if delta == 0 else (1 if delta > 0 else -1)

        consistent = 0
        total = 0
        for dropped in runs:
            sub_runs = [r for r in runs if r != dropped]
            sub_m = [run_values[r] for r in sub_runs if groups[r] == "M"]
            sub_s = [run_values[r] for r in sub_runs if groups[r] == "S"]
            sub_delta = mean(sub_s) - mean(sub_m)
            sub_sign = 0 if sub_delta == 0 else (1 if sub_delta > 0 else -1)
            is_consistent = sub_sign == full_sign
            consistent += int(is_consistent)
            total += 1
            loo_rows.append(
                {
                    "module": module,
                    "dropped_run": dropped,
                    "dropped_group": groups[dropped],
                    "delta_S_minus_M": f"{sub_delta:.12g}",
                    "direction_consistent_with_full": str(is_consistent),
                }
            )

        stat_rows.append(
            {
                "module": module,
                "pathways_defined": str(len(modules[module])),
                "pathways_present": str(module_sizes[module]),
                "n_M": str(len(m_vals)),
                "n_S": str(len(s_vals)),
                "mean_logratio_M": f"{mean(m_vals):.12g}",
                "mean_logratio_S": f"{mean(s_vals):.12g}",
                "median_logratio_M": f"{median(m_vals):.12g}",
                "median_logratio_S": f"{median(s_vals):.12g}",
                "delta_mean_logratio_S_minus_M": f"{delta:.12g}",
                "cliffs_delta_S_vs_M": f"{cliffs_delta(m_vals, s_vals):.12g}",
                "p_exact": f"{p_exact:.12g}",
                "bootstrap95_ci_low": f"{ci_lo:.12g}",
                "bootstrap95_ci_high": f"{ci_hi:.12g}",
                "loo_direction_consistent_n": str(consistent),
                "loo_total": str(total),
            }
        )
    return stat_rows, loo_rows


def pathway_metrics(runs: list[str], rel: dict[str, dict[str, float]]) -> dict[str, dict[str, float]]:
    metrics: dict[str, dict[str, float]] = {}
    for feature, row in rel.items():
        vals = [row.get(run, 0.0) for run in runs]
        m = mean(vals)
        var = mean([(v - m) ** 2 for v in vals]) if vals else 0.0
        prevalence = sum(1 for v in vals if v > 0.0) / len(vals) if vals else 0.0
        zero_fraction = sum(1 for v in vals if v == 0.0) / len(vals) if vals else 0.0
        metrics[feature] = {"mean": m, "var": var, "prevalence": prevalence, "zero_fraction": zero_fraction}
    return metrics


def quantile_bins(metrics: dict[str, dict[str, float]], key: str, n_bins: int) -> dict[str, int]:
    ordered = sorted(metrics, key=lambda f: (metrics[f][key], f))
    bins: dict[str, int] = {}
    if not ordered:
        return bins
    for rank, feature in enumerate(ordered):
        bins[feature] = min(n_bins - 1, int(rank * n_bins / len(ordered)))
    return bins


def random_module_rows(
    runs: list[str],
    rel: dict[str, dict[str, float]],
    modules: dict[str, set[str]],
    groups: dict[str, str],
    pseudocount: float,
    n_random: int,
    seed: int,
) -> list[dict[str, str]]:
    target_features = resolve_module_features(modules.get("overall_fermentation", set()), set(rel))
    target = set(target_features)
    if not target_features:
        return []
    metrics = pathway_metrics(runs, rel)
    mean_bins = quantile_bins(metrics, "mean", 5)
    prev_bins = quantile_bins(metrics, "prevalence", 5)
    var_bins = quantile_bins(metrics, "var", 5)
    zero_bins = quantile_bins(metrics, "zero_fraction", 5)
    by_bin: dict[tuple[int, int, int, int], list[str]] = defaultdict(list)
    by_mean_prev: dict[tuple[int, int], list[str]] = defaultdict(list)
    for feature in sorted(metrics):
        if feature in target:
            continue
        by_bin[(mean_bins[feature], prev_bins[feature], var_bins[feature], zero_bins[feature])].append(feature)
        by_mean_prev[(mean_bins[feature], prev_bins[feature])].append(feature)

    rng = random.Random(seed)
    target_bins = [(mean_bins[f], prev_bins[f], var_bins[f], zero_bins[f]) for f in target_features]
    full_values = score_one_module(set(target), runs, rel, pseudocount)
    full_m = [full_values[r] for r in runs if groups[r] == "M"]
    full_s = [full_values[r] for r in runs if groups[r] == "S"]
    full_delta = mean(full_s) - mean(full_m)
    rows: list[dict[str, str]] = []

    all_pool = [f for f in sorted(metrics) if f not in target]
    for i in range(1, n_random + 1):
        chosen: list[str] = []
        used: set[str] = set()
        exact_bin_matches = 0
        relaxed_mean_prev_matches = 0
        global_fallback_matches = 0
        for bin_key in target_bins:
            candidates = [f for f in by_bin.get(bin_key, []) if f not in used]
            fallback_type = "exact"
            if not candidates:
                mean_prev_key = (bin_key[0], bin_key[1])
                candidates = [f for f in by_mean_prev.get(mean_prev_key, []) if f not in used]
                fallback_type = "relaxed_mean_prevalence"
            if not candidates:
                candidates = [f for f in all_pool if f not in used]
                fallback_type = "global"
            pick = rng.choice(candidates)
            chosen.append(pick)
            used.add(pick)
            if fallback_type == "exact":
                exact_bin_matches += 1
            elif fallback_type == "relaxed_mean_prevalence":
                relaxed_mean_prev_matches += 1
            else:
                global_fallback_matches += 1
        values = score_one_module(set(chosen), runs, rel, pseudocount)
        m_vals = [values[r] for r in runs if groups[r] == "M"]
        s_vals = [values[r] for r in runs if groups[r] == "S"]
        delta = mean(s_vals) - mean(m_vals)
        rows.append(
            {
                "random_module_id": f"random_{i:05d}",
                "n_pathways": str(len(chosen)),
                "delta_mean_logratio_S_minus_M": f"{delta:.12g}",
                "abs_delta_ge_true": str(abs(delta) >= abs(full_delta)),
                "same_direction_as_true": str((delta < 0 and full_delta < 0) or (delta > 0 and full_delta > 0)),
                "p_exact": "not_computed_for_random_specificity",
                "matching_method": "pathway_count_plus_mean_abundance_prevalence_variance_zero_fraction_quintiles",
                "exact_bin_matches": str(exact_bin_matches),
                "relaxed_mean_prevalence_matches": str(relaxed_mean_prev_matches),
                "global_fallback_matches": str(global_fallback_matches),
                "pathways": ";".join(chosen),
            }
        )
    return rows


def quantile(sorted_values: list[float], q: float) -> float:
    if not sorted_values:
        return float("nan")
    if len(sorted_values) == 1:
        return sorted_values[0]
    pos = q * (len(sorted_values) - 1)
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return sorted_values[lo]
    return sorted_values[lo] * (hi - pos) + sorted_values[hi] * (pos - lo)


def empirical_summary(random_rows: list[dict[str, str]], true_stats: list[dict[str, str]]) -> list[dict[str, str]]:
    true = next((r for r in true_stats if r["module"] == "overall_fermentation"), None)
    if not true or not random_rows:
        return []
    true_delta = float(true["delta_mean_logratio_S_minus_M"])
    ge = sum(1 for r in random_rows if r["abs_delta_ge_true"] == "True")
    same = sum(1 for r in random_rows if r["same_direction_as_true"] == "True")
    deltas = [float(r["delta_mean_logratio_S_minus_M"]) for r in random_rows]
    if true_delta < 0:
        directional_extreme = sum(1 for delta in deltas if delta <= true_delta)
    elif true_delta > 0:
        directional_extreme = sum(1 for delta in deltas if delta >= true_delta)
    else:
        directional_extreme = sum(1 for delta in deltas if delta == 0)
    n = len(random_rows)
    return [
        {
            "module": "overall_fermentation",
            "true_delta_mean_logratio_S_minus_M": f"{true_delta:.12g}",
            "n_random_modules": str(n),
            "empirical_p_abs_delta": f"{(ge + 1) / (n + 1):.12g}",
            "empirical_p_directional_delta": f"{(directional_extreme + 1) / (n + 1):.12g}",
            "same_direction_fraction": f"{same / n:.12g}",
            "matching_method": "pathway_count_plus_mean_abundance_prevalence_variance_zero_fraction_quintiles",
        }
    ]


def random_distribution_summary(random_rows: list[dict[str, str]], true_stats: list[dict[str, str]]) -> list[dict[str, str]]:
    true = next((r for r in true_stats if r["module"] == "overall_fermentation"), None)
    if not true or not random_rows:
        return []
    deltas = sorted(float(r["delta_mean_logratio_S_minus_M"]) for r in random_rows)
    abs_deltas = sorted(abs(x) for x in deltas)
    exact = sum(int(r.get("exact_bin_matches", "0")) for r in random_rows)
    relaxed = sum(int(r.get("relaxed_mean_prevalence_matches", "0")) for r in random_rows)
    global_fallback = sum(int(r.get("global_fallback_matches", "0")) for r in random_rows)
    total_slots = exact + relaxed + global_fallback
    return [
        {
            "module": "overall_fermentation",
            "true_delta_mean_logratio_S_minus_M": true["delta_mean_logratio_S_minus_M"],
            "n_random_modules": str(len(random_rows)),
            "random_delta_min": f"{deltas[0]:.12g}",
            "random_delta_q025": f"{quantile(deltas, 0.025):.12g}",
            "random_delta_q25": f"{quantile(deltas, 0.25):.12g}",
            "random_delta_median": f"{quantile(deltas, 0.5):.12g}",
            "random_delta_q75": f"{quantile(deltas, 0.75):.12g}",
            "random_delta_q975": f"{quantile(deltas, 0.975):.12g}",
            "random_delta_max": f"{deltas[-1]:.12g}",
            "random_abs_delta_q95": f"{quantile(abs_deltas, 0.95):.12g}",
            "exact_bin_match_fraction": f"{exact / total_slots:.12g}" if total_slots else "nan",
            "relaxed_mean_prevalence_match_fraction": f"{relaxed / total_slots:.12g}" if total_slots else "nan",
            "global_fallback_fraction": f"{global_fallback / total_slots:.12g}" if total_slots else "nan",
        }
    ]


def comparison_module_rows(stat_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    controls = {
        "bile_acids": "non-SCFA bile-acid transformation comparator",
        "lps_lipidA": "non-fermentation cell-envelope/endotoxin comparator",
        "tryptophan_indole": "non-fermentation amino-acid comparator",
    }
    rows: list[dict[str, str]] = []
    for row in stat_rows:
        if row["module"] not in controls:
            continue
        out = dict(row)
        out["comparison_module_rationale"] = controls[row["module"]]
        rows.append(out)
    return rows


def comparison_module_decisions(stat_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    primary = next((row for row in stat_rows if row["module"] == "overall_fermentation"), None)
    controls = comparison_module_rows(stat_rows)
    if not primary or not controls:
        return []

    primary_delta = float(primary["delta_mean_logratio_S_minus_M"])
    primary_sign = 0 if primary_delta == 0 else (1 if primary_delta > 0 else -1)
    primary_abs = abs(primary_delta)
    primary_p = float(primary["p_exact"])
    primary_loo = int(primary["loo_direction_consistent_n"]) / int(primary["loo_total"])

    rows: list[dict[str, str]] = []
    for control in controls:
        delta = float(control["delta_mean_logratio_S_minus_M"])
        sign = 0 if delta == 0 else (1 if delta > 0 else -1)
        abs_delta = abs(delta)
        p_exact = float(control["p_exact"])
        loo = int(control["loo_direction_consistent_n"]) / int(control["loo_total"])
        same_direction = sign == primary_sign
        equal_or_larger_abs = abs_delta >= primary_abs
        comparable_or_smaller_p = p_exact <= primary_p
        comparable_or_stronger_loo = loo >= primary_loo
        concerning = same_direction and equal_or_larger_abs and comparable_or_smaller_p and comparable_or_stronger_loo
        if concerning:
            interpretation = "stronger_than_primary_on_all_descriptive_checks"
        elif same_direction and equal_or_larger_abs:
            interpretation = "large_but_not_statistically_or_stability_comparable"
        elif same_direction:
            interpretation = "same_direction_weaker_than_primary"
        else:
            interpretation = "opposite_or_unstable_direction"
        rows.append(
            {
                "module": control["module"],
                "comparison_module_rationale": control["comparison_module_rationale"],
                "delta_mean_logratio_S_minus_M": control["delta_mean_logratio_S_minus_M"],
                "p_exact": control["p_exact"],
                "loo_fraction": f"{loo:.12g}",
                "same_direction_as_primary": str(same_direction),
                "equal_or_larger_abs_delta_than_primary": str(equal_or_larger_abs),
                "comparable_or_smaller_p_than_primary": str(comparable_or_smaller_p),
                "comparable_or_stronger_loo_than_primary": str(comparable_or_stronger_loo),
                "compound_concern_flag": str(concerning),
                "interpretation": interpretation,
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise SystemExit(f"No rows to write for {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="PRJDB36442 finalised module reanalysis.")
    parser.add_argument(
        "--merged-pathabundance",
        type=Path,
        default=Path("results/processed/PRJDB36442_humann/merged_pathabundance.tsv.gz"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("metadata/cohorts/PRJDB36442_manifest.tsv"),
    )
    parser.add_argument(
        "--module-definitions",
        type=Path,
        default=Path("analysis_plan/module_definitions.tsv"),
    )
    parser.add_argument("--membership", choices=["conservative", "expanded"], default="conservative")
    parser.add_argument("--pseudocount", type=float, default=1e-9)
    parser.add_argument("--n-boot", type=int, default=5000)
    parser.add_argument("--n-random", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260714)
    parser.add_argument("--out-dir", type=Path, default=Path("results/restructure/PRJDB36442"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = read_tsv(args.manifest)
    groups = {row["run_accession"]: row["group"] for row in manifest}
    modules = read_modules(args.module_definitions, args.membership)
    runs, rel = read_pathway_abundance(args.merged_pathabundance)
    overlapping_runs = [run for run in runs if groups.get(run) in {"M", "S"}]
    if len(overlapping_runs) != 20:
        raise SystemExit(f"Expected 20 overlapping PRJDB36442 runs, found {len(overlapping_runs)}")
    score_rows, values_by_module, present_counts = build_scores(
        overlapping_runs, rel, modules, groups, args.pseudocount
    )
    stat_rows, loo_rows = module_stats(values_by_module, groups, present_counts, modules, args.seed, args.n_boot)
    prefix = args.out_dir / args.membership
    write_tsv(prefix.with_name(f"{args.membership}_module_scores.tsv"), score_rows)
    write_tsv(prefix.with_name(f"{args.membership}_module_stats.tsv"), stat_rows)
    write_tsv(prefix.with_name(f"{args.membership}_leave_one_out.tsv"), loo_rows)
    random_rows = random_module_rows(
        overlapping_runs, rel, modules, groups, args.pseudocount, args.n_random, args.seed
    )
    write_tsv(prefix.with_name(f"{args.membership}_random_modules.tsv"), random_rows)
    summary_rows = empirical_summary(random_rows, stat_rows)
    write_tsv(prefix.with_name(f"{args.membership}_random_module_empirical.tsv"), summary_rows)
    distribution_rows = random_distribution_summary(random_rows, stat_rows)
    write_tsv(prefix.with_name(f"{args.membership}_random_module_distribution.tsv"), distribution_rows)
    control_rows = comparison_module_rows(stat_rows)
    write_tsv(prefix.with_name(f"{args.membership}_mechanistic_comparison_modules.tsv"), control_rows)
    control_decision_rows = comparison_module_decisions(stat_rows)
    write_tsv(prefix.with_name(f"{args.membership}_mechanistic_comparison_summary.tsv"), control_decision_rows)
    print(f"[ok] wrote {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
