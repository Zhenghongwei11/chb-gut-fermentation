#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import random
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def pathway_id(feature: str) -> str:
    return str(feature).split(": ", 1)[0]


def load_modules(path: Path, membership: str) -> dict[str, list[dict[str, str]]]:
    definitions = pd.read_csv(path, sep="\t")
    definitions = definitions[
        (definitions["membership"] == membership)
        & (definitions["lock_status"].isin(["finalised", "draft_seed_not_final", "locked", "draft_seed_not_final_locked"]))
    ].copy()

    modules: dict[str, list[dict[str, str]]] = {}
    for module, df in definitions.groupby("module", sort=True):
        modules[module] = df.to_dict(orient="records")

    fermentation_rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for module, rows in modules.items():
        if not module.startswith("scfa_"):
            continue
        for row in rows:
            if row["pathway_id"] in seen:
                continue
            seen.add(row["pathway_id"])
            fermentation_rows.append(dict(row, module="overall_fermentation"))
    if fermentation_rows:
        modules["overall_fermentation"] = fermentation_rows
    return modules


def prepare_abundance(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t").rename(columns={"# Pathway": "feature"})
    if "feature" not in df.columns:
        raise ValueError("Pathway abundance table must contain feature or # Pathway column")
    df = df[~df["feature"].isin(["UNMAPPED", "UNINTEGRATED"])].copy()
    df = df[~df["feature"].astype(str).str.contains(r"\|", regex=True)].copy()
    df["pathway_id"] = df["feature"].map(pathway_id)
    sample_cols = [c for c in df.columns if c not in {"feature", "pathway_id"}]
    df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    totals = df[sample_cols].sum(axis=0)
    rel = df.copy()
    rel[sample_cols] = rel[sample_cols].div(totals.where(totals > 0, np.nan), axis=1).fillna(0.0)
    return rel.set_index("pathway_id", drop=False)


def cliffs_delta(reference: np.ndarray, test: np.ndarray) -> float:
    if len(reference) == 0 or len(test) == 0:
        return float("nan")
    gt = 0
    lt = 0
    for y in test:
        gt += int(np.sum(y > reference))
        lt += int(np.sum(y < reference))
    return float((gt - lt) / (len(reference) * len(test)))


def bootstrap_ci(reference: np.ndarray, test: np.ndarray, n_boot: int, seed: int) -> tuple[float, float]:
    if len(reference) == 0 or len(test) == 0 or n_boot <= 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    deltas = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        x = rng.choice(reference, size=len(reference), replace=True)
        y = rng.choice(test, size=len(test), replace=True)
        deltas[i] = float(np.mean(y) - np.mean(x))
    return float(np.quantile(deltas, 0.025)), float(np.quantile(deltas, 0.975))


def permutation_p(
    reference: np.ndarray,
    test: np.ndarray,
    n_permutations: int,
    seed: int,
    exact_limit: int = 100000,
) -> tuple[float, str, int]:
    if len(reference) == 0 or len(test) == 0:
        return float("nan"), "not_available", 0
    observed = abs(float(np.mean(test) - np.mean(reference)))
    values = np.concatenate([reference, test])
    n_ref = len(reference)
    n_total = len(values)
    n_combinations = math.comb(n_total, n_ref)
    extreme = 0
    total = 0
    if n_combinations <= exact_limit:
        all_idx = set(range(n_total))
        for ref_idx_tuple in combinations(range(n_total), n_ref):
            ref_idx = set(ref_idx_tuple)
            test_idx = list(all_idx - ref_idx)
            delta = abs(float(np.mean(values[test_idx]) - np.mean(values[list(ref_idx)])))
            extreme += int(delta >= observed - 1e-15)
            total += 1
        return float(extreme / total), "exact", total

    rng = np.random.default_rng(seed)
    for _ in range(n_permutations):
        perm = rng.permutation(n_total)
        ref_idx = perm[:n_ref]
        test_idx = perm[n_ref:]
        delta = abs(float(np.mean(values[test_idx]) - np.mean(values[ref_idx])))
        extreme += int(delta >= observed - 1e-15)
        total += 1
    return float((extreme + 1) / (total + 1)), "monte_carlo", total


def score_pathway_set(rel: pd.DataFrame, sample_ids: list[str], pathway_ids: list[str], pseudocount: float) -> pd.Series:
    module_sum = rel.loc[pathway_ids, sample_ids].sum(axis=0)
    background = (1.0 - module_sum).clip(lower=0.0)
    return np.log(module_sum + pseudocount) - np.log(background + pseudocount)


def score_modules(
    rel: pd.DataFrame,
    metadata: pd.DataFrame,
    modules: dict[str, list[dict[str, str]]],
    membership: str,
    reference_group: str,
    test_group: str,
    contrast_name: str,
    pseudocount: float,
    n_boot: int,
    n_permutations: int,
    seed: int,
    cohort: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample_ids = metadata["sample_id"].tolist()
    score_rows: list[dict[str, object]] = []
    stats_rows: list[dict[str, object]] = []
    member_rows: list[dict[str, object]] = []
    loo_rows: list[dict[str, object]] = []
    map_rows: list[dict[str, object]] = []

    for module, rows in modules.items():
        defined_ids = [row["pathway_id"] for row in rows]
        present_ids = [pid for pid in defined_ids if pid in rel.index]
        for row in rows:
            pid = row["pathway_id"]
            map_rows.append(
                {
                    "cohort": cohort,
                    "membership": membership,
                    "module": module,
                    "pathway_id": pid,
                    "pathway_name": row["pathway_name"],
                    "expected_direction_with_severity": row["expected_direction_with_severity"],
                    "present_in_input": pid in rel.index,
                    "feature": rel.loc[pid, "feature"] if pid in rel.index else "",
                }
            )
        if not present_ids:
            continue

        scores = score_pathway_set(rel, sample_ids, present_ids, pseudocount)
        score_meta = metadata[["sample_id", "group"]].copy()
        score_meta["module_score"] = score_meta["sample_id"].map(scores.to_dict()).astype(float)
        score_meta["relative_module_sum"] = score_meta["sample_id"].map(
            rel.loc[present_ids, sample_ids].sum(axis=0).to_dict()
        ).astype(float)

        for _, row in score_meta.iterrows():
            score_rows.append(
                {
                    "cohort": cohort,
                    "sample_id": row["sample_id"],
                    "group": row["group"],
                    "membership": membership,
                    "module": module,
                    "module_score": row["module_score"],
                    "relative_module_sum": row["relative_module_sum"],
                    "n_pathways_defined": len(defined_ids),
                    "n_pathways_present": len(present_ids),
                    "evidence_tier": "processed_public_matrix_directional_support",
                }
            )

        ref = score_meta.loc[score_meta["group"] == reference_group, "module_score"].to_numpy(dtype=float)
        test = score_meta.loc[score_meta["group"] == test_group, "module_score"].to_numpy(dtype=float)
        if len(ref) and len(test):
            delta = float(np.mean(test) - np.mean(ref))
            p_wilcox = float(mannwhitneyu(ref, test, alternative="two-sided").pvalue)
            p_perm, perm_method, perm_n = permutation_p(ref, test, n_permutations, seed + len(stats_rows))
            ci_low, ci_high = bootstrap_ci(ref, test, n_boot, seed + 1000 + len(stats_rows))
            cliff = cliffs_delta(ref, test)
        else:
            delta = p_wilcox = p_perm = ci_low = ci_high = cliff = float("nan")
            perm_method = "not_available"
            perm_n = 0

        stats_rows.append(
            {
                "cohort": cohort,
                "membership": membership,
                "module": module,
                "contrast": contrast_name,
                "reference_group": reference_group,
                "test_group": test_group,
                "n_reference": len(ref),
                "n_test": len(test),
                "mean_reference": float(np.mean(ref)) if len(ref) else np.nan,
                "mean_test": float(np.mean(test)) if len(test) else np.nan,
                "median_reference": float(np.median(ref)) if len(ref) else np.nan,
                "median_test": float(np.median(test)) if len(test) else np.nan,
                "delta_mean": delta,
                "cliffs_delta_test_vs_reference": cliff,
                "wilcox_p": p_wilcox,
                "permutation_p": p_perm,
                "permutation_method": perm_method,
                "n_permutations": perm_n,
                "bootstrap95_ci_low": ci_low,
                "bootstrap95_ci_high": ci_high,
                "n_pathways_defined": len(defined_ids),
                "n_pathways_present": len(present_ids),
                "evidence_tier": "processed_public_matrix_directional_support",
            }
        )

        for row in rows:
            pid = row["pathway_id"]
            base = {
                "cohort": cohort,
                "membership": membership,
                "module": module,
                "pathway_id": pid,
                "pathway_name": row["pathway_name"],
                "expected_direction_with_severity": row["expected_direction_with_severity"],
                "present_in_input": pid in rel.index,
                "contrast": contrast_name,
                "reference_group": reference_group,
                "test_group": test_group,
            }
            if pid not in rel.index:
                member_rows.append(base)
                continue
            vals = rel.loc[pid, sample_ids]
            tmp = metadata[["sample_id", "group"]].copy()
            tmp["relative_abundance"] = tmp["sample_id"].map(vals.to_dict()).astype(float)
            ref_vals = tmp.loc[tmp["group"] == reference_group, "relative_abundance"].to_numpy(dtype=float)
            test_vals = tmp.loc[tmp["group"] == test_group, "relative_abundance"].to_numpy(dtype=float)
            member_delta = float(np.mean(test_vals) - np.mean(ref_vals))
            expected = row["expected_direction_with_severity"]
            if expected == "decrease":
                expected_met = member_delta < 0
            elif expected == "increase":
                expected_met = member_delta > 0
            else:
                expected_met = np.nan
            member_rows.append(
                {
                    **base,
                    "mean_reference": float(np.mean(ref_vals)),
                    "mean_test": float(np.mean(test_vals)),
                    "delta_mean": member_delta,
                    "expected_direction_met": expected_met,
                    "prevalence": float((vals > 0).mean()),
                }
            )

        full_delta = delta
        for dropped in present_ids:
            kept = [pid for pid in present_ids if pid != dropped]
            if not kept:
                loo_rows.append(
                    {
                        "cohort": cohort,
                        "membership": membership,
                        "module": module,
                        "contrast": contrast_name,
                        "dropped_pathway_id": dropped,
                        "n_pathways_remaining": 0,
                        "delta_mean": np.nan,
                        "same_direction_as_full": np.nan,
                    }
                )
                continue
            loo_score = score_pathway_set(rel, sample_ids, kept, pseudocount)
            ref_loo = loo_score.loc[metadata.loc[metadata["group"] == reference_group, "sample_id"]].to_numpy(dtype=float)
            test_loo = loo_score.loc[metadata.loc[metadata["group"] == test_group, "sample_id"]].to_numpy(dtype=float)
            loo_delta = float(np.mean(test_loo) - np.mean(ref_loo))
            same = (loo_delta < 0 and full_delta < 0) or (loo_delta > 0 and full_delta > 0) or (
                loo_delta == 0 and full_delta == 0
            )
            loo_rows.append(
                {
                    "cohort": cohort,
                    "membership": membership,
                    "module": module,
                    "contrast": contrast_name,
                    "dropped_pathway_id": dropped,
                    "n_pathways_remaining": len(kept),
                    "delta_mean": loo_delta,
                    "same_direction_as_full": same,
                }
            )

    return (
        pd.DataFrame(score_rows),
        pd.DataFrame(stats_rows),
        pd.DataFrame(member_rows),
        pd.DataFrame(loo_rows),
        pd.DataFrame(map_rows),
    )


def quantile_bins(metrics: dict[str, dict[str, float]], key: str, n_bins: int) -> dict[str, int]:
    ordered = sorted(metrics, key=lambda feature: (metrics[feature][key], feature))
    return {feature: min(n_bins - 1, int(rank * n_bins / len(ordered))) for rank, feature in enumerate(ordered)}


def random_module_calibration(
    rel: pd.DataFrame,
    metadata: pd.DataFrame,
    modules: dict[str, list[dict[str, str]]],
    membership: str,
    reference_group: str,
    test_group: str,
    contrast_name: str,
    n_random: int,
    seed: int,
    pseudocount: float,
    cohort: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    target_ids: list[str] = []
    seen: set[str] = set()
    for row in modules.get("overall_fermentation", []):
        pid = row["pathway_id"]
        if pid in rel.index and pid not in seen:
            target_ids.append(pid)
            seen.add(pid)
    if not target_ids:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    sample_ids = metadata["sample_id"].tolist()
    values = rel[sample_ids]
    metrics: dict[str, dict[str, float]] = {}
    for pathway in rel.index.tolist():
        vals = values.loc[pathway].to_numpy(dtype=float)
        mean_value = float(np.mean(vals))
        metrics[pathway] = {
            "mean": mean_value,
            "var": float(np.mean((vals - mean_value) ** 2)),
            "prevalence": float(np.mean(vals > 0.0)),
            "zero_fraction": float(np.mean(vals == 0.0)),
        }

    mean_bins = quantile_bins(metrics, "mean", 5)
    prev_bins = quantile_bins(metrics, "prevalence", 5)
    var_bins = quantile_bins(metrics, "var", 5)
    zero_bins = quantile_bins(metrics, "zero_fraction", 5)
    target = set(target_ids)
    by_bin: dict[tuple[int, int, int, int], list[str]] = {}
    by_mean_prev: dict[tuple[int, int], list[str]] = {}
    for pathway in sorted(metrics):
        if pathway in target:
            continue
        full_key = (mean_bins[pathway], prev_bins[pathway], var_bins[pathway], zero_bins[pathway])
        mean_prev_key = (mean_bins[pathway], prev_bins[pathway])
        by_bin.setdefault(full_key, []).append(pathway)
        by_mean_prev.setdefault(mean_prev_key, []).append(pathway)

    all_pool = [pathway for pathway in sorted(metrics) if pathway not in target]
    if len(all_pool) < len(target_ids):
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    reference = metadata.loc[metadata["group"] == reference_group, "sample_id"].tolist()
    test = metadata.loc[metadata["group"] == test_group, "sample_id"].tolist()
    target_scores = score_pathway_set(rel, sample_ids, target_ids, pseudocount)
    true_delta = float(target_scores.loc[test].mean() - target_scores.loc[reference].mean())

    rng = random.Random(seed)
    target_bins = [(mean_bins[p], prev_bins[p], var_bins[p], zero_bins[p]) for p in target_ids]
    rows: list[dict[str, object]] = []
    for i in range(1, n_random + 1):
        chosen: list[str] = []
        used: set[str] = set()
        exact_matches = 0
        relaxed_matches = 0
        global_matches = 0
        for bin_key in target_bins:
            candidates = [p for p in by_bin.get(bin_key, []) if p not in used]
            fallback = "exact"
            if not candidates:
                candidates = [p for p in by_mean_prev.get((bin_key[0], bin_key[1]), []) if p not in used]
                fallback = "relaxed"
            if not candidates:
                candidates = [p for p in all_pool if p not in used]
                fallback = "global"
            pick = rng.choice(candidates)
            chosen.append(pick)
            used.add(pick)
            exact_matches += int(fallback == "exact")
            relaxed_matches += int(fallback == "relaxed")
            global_matches += int(fallback == "global")

        random_scores = score_pathway_set(rel, sample_ids, chosen, pseudocount)
        delta = float(random_scores.loc[test].mean() - random_scores.loc[reference].mean())
        rows.append(
            {
                "cohort": cohort,
                "membership": membership,
                "random_module_id": f"random_{i:05d}",
                "n_pathways": len(chosen),
                "contrast": contrast_name,
                "delta_mean": delta,
                "abs_delta_ge_true": abs(delta) >= abs(true_delta),
                "same_direction_as_true": (delta < 0 and true_delta < 0) or (delta > 0 and true_delta > 0),
                "matching_method": "pathway_count_plus_mean_abundance_prevalence_variance_zero_fraction_quintiles",
                "exact_bin_matches": exact_matches,
                "relaxed_mean_prevalence_matches": relaxed_matches,
                "global_fallback_matches": global_matches,
                "pathways": ";".join(chosen),
            }
        )

    random_df = pd.DataFrame(rows)
    deltas = random_df["delta_mean"].astype(float).to_numpy()
    n = len(deltas)
    if true_delta < 0:
        directional_extreme = int(np.sum(deltas <= true_delta))
    elif true_delta > 0:
        directional_extreme = int(np.sum(deltas >= true_delta))
    else:
        directional_extreme = int(np.sum(deltas == 0.0))
    total_slots = (
        random_df["exact_bin_matches"].sum()
        + random_df["relaxed_mean_prevalence_matches"].sum()
        + random_df["global_fallback_matches"].sum()
    )
    summary = pd.DataFrame(
        [
            {
                "cohort": cohort,
                "membership": membership,
                "module": "overall_fermentation",
                "contrast": contrast_name,
                "true_delta_mean": true_delta,
                "n_random_modules": n,
                "empirical_p_abs_delta": (int(np.sum(np.abs(deltas) >= abs(true_delta))) + 1) / (n + 1),
                "empirical_p_directional_delta": (directional_extreme + 1) / (n + 1),
                "same_direction_fraction": float(np.mean(random_df["same_direction_as_true"].astype(bool))),
                "matching_method": "pathway_count_plus_mean_abundance_prevalence_variance_zero_fraction_quintiles",
            }
        ]
    )
    distribution = pd.DataFrame(
        [
            {
                "cohort": cohort,
                "membership": membership,
                "module": "overall_fermentation",
                "contrast": contrast_name,
                "true_delta_mean": true_delta,
                "n_random_modules": n,
                "random_delta_min": float(np.min(deltas)),
                "random_delta_q025": float(np.quantile(deltas, 0.025)),
                "random_delta_q25": float(np.quantile(deltas, 0.25)),
                "random_delta_median": float(np.quantile(deltas, 0.5)),
                "random_delta_q75": float(np.quantile(deltas, 0.75)),
                "random_delta_q975": float(np.quantile(deltas, 0.975)),
                "random_delta_max": float(np.max(deltas)),
                "random_abs_delta_q95": float(np.quantile(np.abs(deltas), 0.95)),
                "exact_bin_match_fraction": float(random_df["exact_bin_matches"].sum() / total_slots) if total_slots else np.nan,
                "relaxed_mean_prevalence_match_fraction": float(random_df["relaxed_mean_prevalence_matches"].sum() / total_slots) if total_slots else np.nan,
                "global_fallback_fraction": float(random_df["global_fallback_matches"].sum() / total_slots) if total_slots else np.nan,
            }
        ]
    )
    return random_df, summary, distribution


def direction_consistency(members: pd.DataFrame) -> pd.DataFrame:
    if members.empty:
        return pd.DataFrame()
    rows: list[dict[str, object]] = []
    usable = members[(members["present_in_input"] == True) & (members["expected_direction_with_severity"] == "decrease")]
    for (cohort, membership, module), df in usable.groupby(["cohort", "membership", "module"], sort=True):
        rows.append(
            {
                "cohort": cohort,
                "membership": membership,
                "module": module,
                "n_expected_decrease_present": int(df.shape[0]),
                "n_decreased": int((df["delta_mean"].astype(float) < 0).sum()),
                "fraction_decreased": float((df["delta_mean"].astype(float) < 0).mean()) if df.shape[0] else np.nan,
            }
        )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Score locked modules in a processed public pathway matrix.")
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--pathway", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--modules", type=Path, default=Path("analysis_plan/module_definitions.tsv"))
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--reference-group", required=True)
    parser.add_argument("--test-group", required=True)
    parser.add_argument("--contrast", required=True)
    parser.add_argument("--pseudocount", type=float, default=1e-9)
    parser.add_argument("--n-random", type=int, default=5000)
    parser.add_argument("--n-boot", type=int, default=5000)
    parser.add_argument("--n-permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260714)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rel = prepare_abundance(args.pathway)
    metadata = pd.read_csv(args.metadata, sep="\t")
    metadata = metadata[metadata["group"].isin([args.reference_group, args.test_group])].copy()
    metadata = metadata[metadata["sample_id"].isin(rel.columns)].copy()
    if metadata.empty:
        raise ValueError("No metadata rows remain after filtering to requested groups and pathway sample columns")

    all_scores = []
    all_stats = []
    all_members = []
    all_loo = []
    all_maps = []
    all_random = []
    all_random_summary = []
    all_random_distribution = []

    for membership in ["conservative", "expanded"]:
        modules = load_modules(args.modules, membership)
        scores, stats, members, loo, mapping = score_modules(
            rel=rel,
            metadata=metadata,
            modules=modules,
            membership=membership,
            reference_group=args.reference_group,
            test_group=args.test_group,
            contrast_name=args.contrast,
            pseudocount=args.pseudocount,
            n_boot=args.n_boot,
            n_permutations=args.n_permutations,
            seed=args.seed + (0 if membership == "conservative" else 100000),
            cohort=args.cohort,
        )
        random_rows, random_summary, random_distribution = random_module_calibration(
            rel=rel,
            metadata=metadata,
            modules=modules,
            membership=membership,
            reference_group=args.reference_group,
            test_group=args.test_group,
            contrast_name=args.contrast,
            n_random=args.n_random,
            seed=args.seed + (0 if membership == "conservative" else 1),
            pseudocount=args.pseudocount,
            cohort=args.cohort,
        )
        all_scores.append(scores)
        all_stats.append(stats)
        all_members.append(members)
        all_loo.append(loo)
        all_maps.append(mapping)
        all_random.append(random_rows)
        all_random_summary.append(random_summary)
        all_random_distribution.append(random_distribution)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    scores = pd.concat(all_scores, ignore_index=True)
    stats = pd.concat(all_stats, ignore_index=True)
    members = pd.concat(all_members, ignore_index=True)
    loo = pd.concat(all_loo, ignore_index=True)
    mapping = pd.concat(all_maps, ignore_index=True)
    random_modules = pd.concat(all_random, ignore_index=True)
    random_summary = pd.concat(all_random_summary, ignore_index=True)
    random_distribution = pd.concat(all_random_distribution, ignore_index=True)
    direction = direction_consistency(members)

    scores.to_csv(args.out_dir / "module_scores.tsv", sep="\t", index=False)
    stats.to_csv(args.out_dir / "module_binary_contrasts.tsv", sep="\t", index=False)
    members.to_csv(args.out_dir / "pathway_member_binary_stats.tsv", sep="\t", index=False)
    loo.to_csv(args.out_dir / "leave_one_pathway_out.tsv", sep="\t", index=False)
    mapping.to_csv(args.out_dir / "locked_module_feature_mapping.tsv", sep="\t", index=False)
    direction.to_csv(args.out_dir / "pathway_direction_consistency.tsv", sep="\t", index=False)
    random_modules.to_csv(args.out_dir / "random_module_distribution.tsv", sep="\t", index=False)
    random_summary.to_csv(args.out_dir / "random_module_empirical.tsv", sep="\t", index=False)
    random_distribution.to_csv(args.out_dir / "random_module_distribution_summary.tsv", sep="\t", index=False)

    print(f"[ok] wrote processed binary module outputs to {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
