#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import itertools
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


def clean_run_column(column: str) -> str:
    column = re.sub(r"_(Abundance|Coverage)(-RPKs)?$", "", column)
    column = re.sub(r"\.(pathabundance|pathcoverage)$", "", column)
    return column


def pathway_id(feature: str) -> str:
    return str(feature).split(": ", 1)[0]


def exact_permutation_p(values: list[float], groups: list[str]) -> float:
    s_count = sum(1 for group in groups if group == "S")
    m_count = len(groups) - s_count
    if not s_count or not m_count:
        return float("nan")
    obs_s = np.mean([value for value, group in zip(values, groups) if group == "S"])
    obs_m = np.mean([value for value, group in zip(values, groups) if group == "M"])
    obs = obs_s - obs_m
    total_sum = sum(values)
    extreme = 0
    total = 0
    for subset in itertools.combinations(range(len(values)), s_count):
        s_sum = sum(values[i] for i in subset)
        diff = (s_sum / s_count) - ((total_sum - s_sum) / m_count)
        if abs(diff) >= abs(obs) - 1e-15:
            extreme += 1
        total += 1
    return extreme / total


def load_abundance(path: Path, runs: list[str]) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        df = pd.read_csv(handle, sep="\t")
    df = df.rename(columns={df.columns[0]: "feature"})
    df = df[~df["feature"].isin(["UNMAPPED", "UNINTEGRATED"])].copy()
    df = df[~df["feature"].str.contains(r"\|", regex=True)].copy()
    df["pathway_id"] = df["feature"].map(pathway_id)
    renamed = {column: clean_run_column(column) for column in df.columns if column not in {"feature", "pathway_id"}}
    df = df.rename(columns=renamed)
    missing = [run for run in runs if run not in df.columns]
    if missing:
        raise SystemExit(f"Abundance table missing runs: {', '.join(missing)}")
    df[runs] = df[runs].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    totals = df[runs].sum(axis=0)
    df[runs] = df[runs].div(totals.where(totals > 0, np.nan), axis=1).fillna(0.0)
    return df.set_index("pathway_id", drop=False)


def load_modules(path: Path, membership: str) -> dict[str, list[dict[str, str]]]:
    definitions = pd.read_csv(path, sep="\t")
    definitions = definitions[
        (definitions["membership"] == membership)
        & (definitions["lock_status"].isin(["finalised", "draft_seed_not_final", "locked", "draft_seed_not_final_locked"]))
    ].copy()
    modules: dict[str, list[dict[str, str]]] = {
        module: df.to_dict(orient="records")
        for module, df in definitions.groupby("module", sort=True)
    }
    seen: set[str] = set()
    fermentation_rows: list[dict[str, str]] = []
    for module, rows in modules.items():
        if not module.startswith("scfa_"):
            continue
        for row in rows:
            pid = row["pathway_id"]
            if pid in seen:
                continue
            seen.add(pid)
            fermentation_rows.append(dict(row, module="overall_fermentation"))
    modules["overall_fermentation"] = fermentation_rows
    return modules


def module_score(rel: pd.DataFrame, runs: list[str], pathway_ids: list[str], pseudocount: float) -> pd.Series:
    present = [pid for pid in pathway_ids if pid in rel.index]
    if not present:
        return pd.Series({run: math.log(pseudocount) - math.log(1.0 + pseudocount) for run in runs})
    module_sum = rel.loc[present, runs].sum(axis=0)
    background = (1.0 - module_sum).clip(lower=0.0)
    return np.log(module_sum + pseudocount) - np.log(background + pseudocount)


def contrast_stats(values: pd.Series, groups: pd.Series) -> dict[str, float]:
    m_vals = values.loc[groups == "M"].astype(float)
    s_vals = values.loc[groups == "S"].astype(float)
    ordered_runs = values.index.tolist()
    return {
        "mean_M": float(m_vals.mean()),
        "mean_S": float(s_vals.mean()),
        "delta_S_minus_M": float(s_vals.mean() - m_vals.mean()),
        "prevalence_overall": float((values > 0).mean()),
        "prevalence_M": float((m_vals > 0).mean()),
        "prevalence_S": float((s_vals > 0).mean()),
        "exact_p": exact_permutation_p(values.loc[ordered_runs].astype(float).tolist(), groups.loc[ordered_runs].tolist()),
    }


def member_stats(rel: pd.DataFrame, metadata: pd.DataFrame, modules: dict[str, list[dict[str, str]]], membership: str) -> pd.DataFrame:
    runs = metadata["run_accession"].tolist()
    groups = metadata.set_index("run_accession")["group"]
    rows: list[dict[str, object]] = []
    module_delta_sums: dict[str, float] = {}
    pending: list[dict[str, object]] = []
    for module, module_rows in modules.items():
        seen: set[str] = set()
        for row in module_rows:
            pid = row["pathway_id"]
            if pid in seen:
                continue
            seen.add(pid)
            base = {
                "membership": membership,
                "module": module,
                "pathway_id": pid,
                "pathway_name": row["pathway_name"],
                "expected_direction_with_severity": row["expected_direction_with_severity"],
                "present_in_input": pid in rel.index,
            }
            if pid not in rel.index:
                pending.append(base)
                continue
            values = rel.loc[pid, runs].astype(float)
            stats = contrast_stats(values, groups)
            expected = row["expected_direction_with_severity"]
            base.update(stats)
            base["expected_direction_met"] = bool(stats["delta_S_minus_M"] < 0) if expected == "decrease" else np.nan
            pending.append(base)
            module_delta_sums[module] = module_delta_sums.get(module, 0.0) + stats["delta_S_minus_M"]
    for row in pending:
        delta = row.get("delta_S_minus_M", np.nan)
        denom = module_delta_sums.get(str(row["module"]), np.nan)
        row["contribution_to_module_delta"] = float(delta / denom) if pd.notna(delta) and denom and not np.isnan(denom) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def leave_one_pathway_out(
    rel: pd.DataFrame,
    metadata: pd.DataFrame,
    modules: dict[str, list[dict[str, str]]],
    membership: str,
    pseudocount: float,
) -> pd.DataFrame:
    runs = metadata["run_accession"].tolist()
    groups = metadata.set_index("run_accession")["group"]
    target = []
    seen: set[str] = set()
    for row in modules["overall_fermentation"]:
        pid = row["pathway_id"]
        if pid in rel.index and pid not in seen:
            target.append(pid)
            seen.add(pid)
    full = module_score(rel, runs, target, pseudocount)
    full_stats = contrast_stats(full, groups)
    rows: list[dict[str, object]] = []
    for dropped in target:
        kept = [pid for pid in target if pid != dropped]
        values = module_score(rel, runs, kept, pseudocount)
        stats = contrast_stats(values, groups)
        rows.append(
            {
                "membership": membership,
                "module": "overall_fermentation",
                "dropped_pathway_id": dropped,
                "n_pathways_remaining": len(kept),
                "full_delta_S_minus_M": full_stats["delta_S_minus_M"],
                "loo_delta_S_minus_M": stats["delta_S_minus_M"],
                "loo_exact_p": stats["exact_p"],
                "direction_consistent_with_full": np.sign(stats["delta_S_minus_M"]) == np.sign(full_stats["delta_S_minus_M"]),
            }
        )
    return pd.DataFrame(rows)


def direction_summary(member_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    present = member_df[member_df["present_in_input"] == True].copy()
    for (membership, module), df in present.groupby(["membership", "module"], sort=True):
        expected_decrease = df[df["expected_direction_with_severity"] == "decrease"]
        rows.append(
            {
                "membership": membership,
                "module": module,
                "n_pathways_present": int(df.shape[0]),
                "n_expected_decrease": int(expected_decrease.shape[0]),
                "n_decreased_S_minus_M": int((expected_decrease["delta_S_minus_M"] < 0).sum()) if not expected_decrease.empty else 0,
                "fraction_expected_decrease_met": float((expected_decrease["delta_S_minus_M"] < 0).mean()) if not expected_decrease.empty else np.nan,
                "module_interpretation": "direction_consistency_summary",
            }
        )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="PRJDB36442 pathway-level sensitivity for finalised modules.")
    parser.add_argument("--merged-pathabundance", type=Path, default=Path("results/processed/PRJDB36442_humann/merged_pathabundance.tsv.gz"))
    parser.add_argument("--manifest", type=Path, default=Path("metadata/cohorts/PRJDB36442_manifest.tsv"))
    parser.add_argument("--module-definitions", type=Path, default=Path("analysis_plan/module_definitions.tsv"))
    parser.add_argument("--out-dir", type=Path, default=Path("results/restructure/PRJDB36442"))
    parser.add_argument("--pseudocount", type=float, default=1e-9)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    metadata = pd.read_csv(args.manifest, sep="\t")
    metadata = metadata[metadata["group"].isin(["M", "S"])].copy()
    runs = metadata["run_accession"].tolist()
    rel = load_abundance(args.merged_pathabundance, runs)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    all_members = []
    all_loo = []
    for membership in ["conservative", "expanded"]:
        modules = load_modules(args.module_definitions, membership)
        members = member_stats(rel, metadata, modules, membership)
        loo = leave_one_pathway_out(rel, metadata, modules, membership, args.pseudocount)
        members.to_csv(args.out_dir / f"{membership}_pathway_member_stats.tsv", sep="\t", index=False)
        loo.to_csv(args.out_dir / f"{membership}_leave_one_pathway_out.tsv", sep="\t", index=False)
        all_members.append(members)
        all_loo.append(loo)

    member_df = pd.concat(all_members, ignore_index=True)
    loo_df = pd.concat(all_loo, ignore_index=True)
    direction = direction_summary(member_df)
    member_df.to_csv(args.out_dir / "pathway_member_stats.tsv", sep="\t", index=False)
    loo_df.to_csv(args.out_dir / "leave_one_pathway_out.tsv", sep="\t", index=False)
    direction.to_csv(args.out_dir / "pathway_direction_consistency.tsv", sep="\t", index=False)
    print(f"[ok] wrote pathway sensitivity outputs to {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
