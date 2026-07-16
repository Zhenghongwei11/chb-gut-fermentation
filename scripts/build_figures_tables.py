#!/usr/bin/env python3
from __future__ import annotations

import gzip
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
FIG = OUT / "figures"
TAB = OUT / "tables"
RES = OUT / "results"
PSEUDOCOUNT = 1e-9

MODULE_LABELS = {
    "overall_fermentation": "Selected fermentation composite",
    "scfa_acetate": "Acetate/lactate overlap",
    "scfa_lactate_succinate": "Acetate/lactate overlap",
    "scfa_butyrate": "Butyrate",
    "bile_acids": "Bile-acid transformation",
    "lps_lipidA": "LPS/lipid A biosynthesis",
    "tryptophan_indole": "Tryptophan biosynthesis",
}


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def fmt_num(value: object, digits: int = 3) -> str:
    value = float(value)
    if math.isnan(value):
        return ""
    return f"{value:.{digits}f}"


def fmt_p(value: object) -> str:
    value = float(value)
    if value < 0.001:
        return f"{value:.2e}"
    return f"{value:.3f}"


def row(df: pd.DataFrame, module: str, membership: str | None = None) -> pd.Series:
    out = df[df["module"] == module]
    if membership is not None and "membership" in out.columns:
        out = out[out["membership"] == membership]
    if out.empty:
        raise KeyError(module)
    return out.iloc[0]


def prj_delta(r: pd.Series) -> float:
    return float(r["delta_mean_logratio_S_minus_M"])


def ext_delta(r: pd.Series) -> float:
    return float(r["delta_mean"])


def prj_p(r: pd.Series) -> float:
    return float(r["p_exact"])


def ext_p(r: pd.Series) -> float:
    return float(r["permutation_p"])


def feature_to_id(feature: str) -> str:
    return str(feature).split("|", 1)[0].split(":", 1)[0]


def draw_panel_label(ax, label: str) -> None:
    ax.text(-0.12, 1.05, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")


def build_figure_1() -> None:
    fig, ax = plt.subplots(figsize=(11.2, 6.2))
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.5, 0.94, "Study cohorts and phenotype definitions", ha="center", fontsize=16, fontweight="bold")
    ax.text(0.5, 0.89, "A discovery-derived pathway composite was finalised before external-cohort evaluation.", ha="center", fontsize=10.5)
    boxes = [
        (0.06, 0.59, 0.26, 0.22, "#245766", "CHB discovery cohort\nTreatment-naive CHB\nM: G0-1 and S0-1, n=9\nS: G>1 and/or S>1, n=11"),
        (0.37, 0.59, 0.26, 0.22, "#2c8f83", "NAFLD fibrosis cohort\nBiopsy-proven NAFLD\nF0-F2, n=72\nF3-F4, n=14"),
        (0.68, 0.59, 0.26, 0.22, "#6b7a8f", "Cirrhosis context cohort\nHealthy controls, n=114\nCirrhosis, n=123"),
        (0.10, 0.27, 0.22, 0.16, "#eef4f7", "Discovery contrast\nSelected composite\nand submodules"),
        (0.39, 0.27, 0.22, 0.16, "#eef4f7", "Unadjusted comparison\nSelected composite\nand sensitivity"),
        (0.68, 0.27, 0.22, 0.16, "#eef4f7", "Context contrast\nSelected composite\nand heterogeneity"),
    ]
    for x, y, w, h, color, text in boxes:
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=color, edgecolor="#263238", linewidth=1.0))
        tcolor = "white" if color.startswith("#2") or color.startswith("#6") else "#1f2933"
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=10.1, linespacing=1.25, color=tcolor)
    for x in [0.19, 0.50, 0.81]:
        ax.annotate("", xy=(x, 0.44), xytext=(x, 0.59), arrowprops=dict(arrowstyle="->", lw=1.2, color="#263238"))
    fig.savefig(FIG / "Figure_1.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Figure_1.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_figure_2() -> None:
    scores = read_tsv(RES / "prjdb36442" / "conservative" / "conservative_module_scores.tsv")
    stats = read_tsv(RES / "prjdb36442" / "conservative" / "conservative_module_stats.tsv")
    plot_scores = scores[scores["module"] == "overall_fermentation"].copy()
    rows = [
        ("Selected fermentation composite", row(stats, "overall_fermentation")),
        ("Acetate/lactate overlap", row(stats, "scfa_acetate")),
        ("Butyrate", row(stats, "scfa_butyrate")),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11.3, 4.7), gridspec_kw={"width_ratios": [1.0, 1.25]})
    vals = [plot_scores.loc[plot_scores["group"] == g, "log_ratio_score"].to_numpy() for g in ["M", "S"]]
    axes[0].boxplot(vals, tick_labels=["M\nn=9", "S\nn=11"], patch_artist=True, boxprops=dict(facecolor="#dbeafe"), medianprops=dict(color="#111827", lw=1.5))
    for i, arr in enumerate(vals, start=1):
        axes[0].scatter([i] * len(arr), arr, s=28, color="#245766", edgecolor="white", zorder=3)
    axes[0].set_title("Selected fermentation composite score")
    axes[0].set_ylabel("log-ratio module score")
    axes[0].axhline(0, color="#cbd5e1", lw=0.8)
    draw_panel_label(axes[0], "A")
    labels = [r[0] for r in rows]
    deltas = [prj_delta(r[1]) for r in rows]
    lows = [float(r[1]["bootstrap95_ci_low"]) for r in rows]
    highs = [float(r[1]["bootstrap95_ci_high"]) for r in rows]
    y = list(range(len(rows)))[::-1]
    axes[1].errorbar(deltas, y, xerr=[[d - l for d, l in zip(deltas, lows)], [h - d for d, h in zip(deltas, highs)]], fmt="o", color="#245766", ecolor="#245766", capsize=4)
    axes[1].axvline(0, color="#9ca3af", lw=1)
    axes[1].set_yticks(y, labels)
    axes[1].set_xlabel("S minus M mean difference (95% CI)")
    axes[1].set_title("Discovery-cohort effect estimates")
    for yy, high, (_, stat) in zip(y, highs, rows):
        axes[1].text(high + 0.05, yy, f"P={fmt_p(prj_p(stat))}", va="center", fontsize=9)
    draw_panel_label(axes[1], "B")
    fig.suptitle("CHB discovery analysis", y=1.02, fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "Figure_2.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Figure_2.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_figure_3() -> None:
    scores = read_tsv(RES / "loombar2017" / "module_scores.tsv")
    stats = read_tsv(RES / "loombar2017" / "module_binary_contrasts.tsv")
    plot_scores = scores[(scores["membership"] == "conservative") & (scores["module"] == "overall_fermentation")]
    rows = [
        ("Selected fermentation composite", row(stats, "overall_fermentation", "conservative")),
        ("Acetate/lactate overlap", row(stats, "scfa_acetate", "conservative")),
        ("Butyrate", row(stats, "scfa_butyrate", "conservative")),
        ("Selected composite\nexpanded definition", row(stats, "overall_fermentation", "expanded")),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.8), gridspec_kw={"width_ratios": [1.1, 1.1]})
    vals = [plot_scores.loc[plot_scores["group"] == g, "module_score"].to_numpy() for g in ["F0_F2", "F3_F4"]]
    axes[0].boxplot(vals, tick_labels=["F0-F2\nn=72", "F3-F4\nn=14"], patch_artist=True, boxprops=dict(facecolor="#d1fae5"), medianprops=dict(color="#111827", lw=1.5))
    for i, arr in enumerate(vals, start=1):
        axes[0].scatter([i] * len(arr), arr, s=20, alpha=0.7, color="#2c8f83", edgecolor="white", linewidth=0.3, zorder=3)
    axes[0].set_ylabel("log-ratio module score")
    axes[0].set_title("Selected fermentation composite")
    draw_panel_label(axes[0], "A")
    labels = [r[0] for r in rows]
    rows_for_plot = [r[1] for r in rows]
    deltas = [ext_delta(r) for r in rows_for_plot]
    lows = [float(r["bootstrap95_ci_low"]) for r in rows_for_plot]
    highs = [float(r["bootstrap95_ci_high"]) for r in rows_for_plot]
    y = list(range(len(labels)))[::-1]
    axes[1].errorbar(deltas, y, xerr=[[d - l for d, l in zip(deltas, lows)], [h - d for d, h in zip(deltas, highs)]], fmt="o", color="#2c8f83", ecolor="#2c8f83", capsize=4)
    axes[1].axvline(0, color="#9ca3af", lw=1)
    axes[1].set_yticks(y, labels)
    axes[1].set_xlabel("F3-F4 minus F0-F2 mean difference (95% CI)")
    axes[1].set_title("Finalised modules and expanded sensitivity")
    for yy, high, stat in zip(y, highs, rows_for_plot):
        axes[1].text(high + 0.04, yy, f"P={fmt_p(ext_p(stat))}", va="center", fontsize=9)
    draw_panel_label(axes[1], "B")
    fig.suptitle("NAFLD biopsy-defined fibrosis comparison and pathway sensitivity", y=1.02, fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "Figure_3.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Figure_3.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_figure_4() -> None:
    prj = row(read_tsv(RES / "prjdb36442" / "conservative" / "conservative_module_stats.tsv"), "overall_fermentation")
    loomba = row(read_tsv(RES / "loombar2017" / "module_binary_contrasts.tsv"), "overall_fermentation", "conservative")
    qin = row(read_tsv(RES / "qinn2014" / "module_binary_contrasts.tsv"), "overall_fermentation", "conservative")
    rows = [
        ("CHB discovery cohort\nS vs M", prj_delta(prj), float(prj["bootstrap95_ci_low"]), float(prj["bootstrap95_ci_high"]), prj_p(prj)),
        ("NAFLD fibrosis cohort\nF3-F4 vs F0-F2", ext_delta(loomba), float(loomba["bootstrap95_ci_low"]), float(loomba["bootstrap95_ci_high"]), ext_p(loomba)),
        ("Cirrhosis context cohort\ncirrhosis vs healthy", ext_delta(qin), float(qin["bootstrap95_ci_low"]), float(qin["bootstrap95_ci_high"]), ext_p(qin)),
    ]
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    y = list(range(len(rows)))[::-1]
    for yy, item, color in zip(y, rows, ["#245766", "#2c8f83", "#6b7a8f"]):
        _, d, low, high, pval = item
        ax.errorbar(d, yy, xerr=[[d - low], [high - d]], fmt="o", color=color, ecolor=color, capsize=4, ms=7)
        ax.text(high + 0.06, yy, f"P={fmt_p(pval)}", va="center", fontsize=9.5)
    ax.axvline(0, color="#9ca3af", lw=1)
    ax.set_yticks(y, [r[0] for r in rows])
    ax.set_xlabel("Test minus reference mean difference in selected composite score (95% CI)")
    ax.set_title("Selected fermentation-composite effects across analysed cohorts", fontsize=14, fontweight="bold")
    ax.text(0.02, -0.28, "Cohort-specific estimates; different pathway universes; not pooled.", transform=ax.transAxes, fontsize=9.2, color="#4b5563")
    fig.tight_layout()
    fig.savefig(FIG / "Figure_4.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Figure_4.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_supplementary_figure_1() -> None:
    stats = read_tsv(RES / "prjdb36442" / "conservative" / "conservative_module_stats.tsv")
    rand = read_tsv(RES / "prjdb36442" / "conservative" / "conservative_random_modules.tsv")
    rand_emp = row(read_tsv(RES / "prjdb36442" / "conservative" / "conservative_random_module_empirical.tsv"), "overall_fermentation")
    modules = ["overall_fermentation", "bile_acids", "lps_lipidA", "tryptophan_indole"]
    rows = [(MODULE_LABELS[m], row(stats, m)) for m in modules]
    fig, axes = plt.subplots(1, 2, figsize=(11.3, 4.7), gridspec_kw={"width_ratios": [1.1, 1.2]})
    deltas = [prj_delta(r[1]) for r in rows]
    lows = [float(r[1]["bootstrap95_ci_low"]) for r in rows]
    highs = [float(r[1]["bootstrap95_ci_high"]) for r in rows]
    y = list(range(len(rows)))[::-1]
    for yy, d, low, high, color in zip(y, deltas, lows, highs, ["#245766", "#7c8798", "#7c8798", "#7c8798"]):
        axes[0].errorbar(d, yy, xerr=[[d - low], [high - d]], fmt="o", color=color, ecolor=color, capsize=4)
    axes[0].axvline(0, color="#9ca3af", lw=1)
    axes[0].set_yticks(y, [r[0] for r in rows])
    axes[0].set_xlabel("S minus M mean difference (95% CI)")
    axes[0].set_title("Exploratory mechanistic comparison modules", fontsize=11)
    draw_panel_label(axes[0], "A")
    col = "delta_mean_logratio_S_minus_M"
    axes[1].hist(rand[col].astype(float), bins=35, color="#d8dee9", edgecolor="white")
    axes[1].axvline(prj_delta(row(stats, "overall_fermentation")), color="#c2410c", lw=2, label="Observed selected composite")
    axes[1].axvline(rand[col].astype(float).median(), color="#334155", lw=1.5, ls="--", label="Random-module median")
    axes[1].set_xlabel("Random-module S minus M difference")
    axes[1].set_ylabel("Count")
    axes[1].set_title(f"Matched random modules\nabsolute empirical P={fmt_p(rand_emp['empirical_p_abs_delta'])}", fontsize=11)
    axes[1].legend(frameon=False, fontsize=9)
    draw_panel_label(axes[1], "B")
    fig.suptitle("Calibration of the CHB discovery selected fermentation-composite result", y=1.02, fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "Supplementary_Figure_1.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Supplementary_Figure_1.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def table_1() -> pd.DataFrame:
    return pd.DataFrame(
        [
            ["CHB discovery", "Treatment-naive CHB", "M: Scheuer G0-1 and S0-1; n=9", "S: G>1 and/or S>1; n=11", "Binary biopsy-defined M/S group", "Sample-level G, S, ALT, AST, BMI, HBV DNA, sex, and medication exposure unavailable from public metadata used here"],
            ["NAFLD fibrosis", "Biopsy-proven NAFLD", "F0-F2; n=72", "F3-F4; n=14", "Fibrosis stage available; age category fully collinear with fibrosis group", "Read count, read bases, platform, extraction kit; sex, BMI, diabetes unavailable"],
            ["Cirrhosis context", "Cirrhosis case-control cohort", "Healthy controls; n=114", "Cirrhosis; n=123", "Disease status available; not a disease-internal histology contrast", "Detailed patient-level cirrhosis severity covariates unavailable in this analysis"],
        ],
        columns=["cohort", "disease_setting", "reference_group", "test_group", "phenotype_data", "covariates_available"],
    )


def table_2() -> pd.DataFrame:
    prj = row(read_tsv(RES / "prjdb36442" / "conservative" / "conservative_module_stats.tsv"), "overall_fermentation")
    loomba = row(read_tsv(RES / "loombar2017" / "module_binary_contrasts.tsv"), "overall_fermentation", "conservative")
    qin = row(read_tsv(RES / "qinn2014" / "module_binary_contrasts.tsv"), "overall_fermentation", "conservative")
    return pd.DataFrame(
        [
            ["CHB discovery", "M", "S", 9, 11, prj_delta(prj), prj["cliffs_delta_S_vs_M"], prj_p(prj), prj["bootstrap95_ci_low"], prj["bootstrap95_ci_high"], "Discovery-stage"],
            ["NAFLD fibrosis", "F0-F2", "F3-F4", 72, 14, ext_delta(loomba), loomba["cliffs_delta_test_vs_reference"], ext_p(loomba), loomba["bootstrap95_ci_low"], loomba["bootstrap95_ci_high"], "Unadjusted histologic comparison"],
            ["Cirrhosis context", "Healthy", "Cirrhosis", 114, 123, ext_delta(qin), qin["cliffs_delta_test_vs_reference"], ext_p(qin), qin["bootstrap95_ci_low"], qin["bootstrap95_ci_high"], "Context"],
        ],
        columns=["cohort", "reference_group", "test_group", "reference_n", "test_n", "delta", "cliffs_delta", "p_value", "ci_low", "ci_high", "interpretation"],
    )


def unmapped_unintegrated() -> pd.DataFrame:
    rows = []
    resources = [
        ("PRJDB36442", ROOT / "data" / "prjdb36442" / "merged_pathabundance.tsv.gz"),
        ("LoombaR_2017", ROOT / "data" / "loombar2017" / "pathway_abundance_unstratified.tsv"),
        ("QinN_2014", ROOT / "data" / "qinn2014" / "pathway_abundance_unstratified.tsv"),
    ]
    for cohort, path in resources:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as handle:
            raw = pd.read_csv(handle, sep="\t")
        feature = raw.iloc[:, 0].astype(str).map(feature_to_id)
        values = raw.iloc[:, 1:].astype(float)
        total = values.sum(axis=0).replace(0, np.nan)
        for label in ["UNMAPPED", "UNINTEGRATED"]:
            mask = feature.eq(label)
            if mask.any():
                frac = values.loc[mask].sum(axis=0) / total
                rows.append([cohort, label, int(mask.sum()), len(frac), frac.median(), frac.quantile(0.25), frac.quantile(0.75), frac.min(), frac.max()])
            else:
                rows.append([cohort, label, 0, values.shape[1], np.nan, np.nan, np.nan, np.nan, np.nan])
    return pd.DataFrame(rows, columns=["cohort", "feature_class", "feature_rows_present", "sample_n", "median_fraction", "iqr_low", "iqr_high", "minimum", "maximum"])


def write_workbook() -> None:
    workbook = TAB / "supplementary_tables.xlsx"
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        table_1().to_excel(writer, sheet_name="Table 1 Cohorts", index=False)
        table_2().to_excel(writer, sheet_name="Table 2 Effects", index=False)
        read_tsv(ROOT / "analysis" / "module_definitions.tsv").to_excel(writer, sheet_name="Module definitions", index=False)
        read_tsv(RES / "prjdb36442" / "conservative" / "conservative_module_stats.tsv").to_excel(writer, sheet_name="PRJDB conservative", index=False)
        read_tsv(RES / "prjdb36442" / "expanded" / "expanded_module_stats.tsv").to_excel(writer, sheet_name="PRJDB expanded", index=False)
        read_tsv(RES / "prjdb36442" / "sensitivity" / "pathway_member_stats.tsv").to_excel(writer, sheet_name="PRJDB pathway members", index=False)
        read_tsv(RES / "prjdb36442" / "sensitivity" / "leave_one_pathway_out.tsv").to_excel(writer, sheet_name="PRJDB leave one pathway", index=False)
        read_tsv(RES / "prjdb36442" / "conservative" / "conservative_random_module_empirical.tsv").to_excel(writer, sheet_name="PRJDB random empirical", index=False)
        read_tsv(RES / "loombar2017" / "module_binary_contrasts.tsv").to_excel(writer, sheet_name="Loomba effects", index=False)
        read_tsv(RES / "loombar2017" / "pathway_member_binary_stats.tsv").to_excel(writer, sheet_name="Loomba pathway members", index=False)
        read_tsv(RES / "loombar2017" / "leave_one_pathway_out.tsv").to_excel(writer, sheet_name="Loomba leave one pathway", index=False)
        read_tsv(RES / "loombar2017" / "random_module_empirical.tsv").to_excel(writer, sheet_name="Loomba random empirical", index=False)
        read_tsv(RES / "qinn2014" / "module_binary_contrasts.tsv").to_excel(writer, sheet_name="Qin effects", index=False)
        read_tsv(RES / "qinn2014" / "pathway_member_binary_stats.tsv").to_excel(writer, sheet_name="Qin pathway members", index=False)
        unmapped_unintegrated().to_excel(writer, sheet_name="Unmapped summary", index=False)
    from openpyxl import load_workbook

    wb = load_workbook(workbook)
    fill = PatternFill("solid", fgColor="1F4E78")
    font = Font(color="FFFFFF", bold=True)
    for ws in wb.worksheets:
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions
        for cell in ws[1]:
            cell.fill = fill
            cell.font = font
            cell.alignment = Alignment(wrap_text=True, vertical="center")
        for idx in range(1, ws.max_column + 1):
            col = get_column_letter(idx)
            width = 10
            for cell in ws[col]:
                value = "" if cell.value is None else str(cell.value)
                width = max(width, min(56, len(value) + 2))
                cell.alignment = Alignment(wrap_text=True, vertical="top")
            ws.column_dimensions[col].width = width
    wb.save(workbook)
    table_1().to_csv(TAB / "Table_1_cohort_characteristics.csv", index=False)
    table_2().to_csv(TAB / "Table_2_selected_composite_effects.csv", index=False)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    build_figure_1()
    build_figure_2()
    build_figure_3()
    build_figure_4()
    build_supplementary_figure_1()
    write_workbook()


if __name__ == "__main__":
    main()
