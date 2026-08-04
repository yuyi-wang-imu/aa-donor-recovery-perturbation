from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Mandatory publication settings: editable SVG text and consistent sans-serif type.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"

import matplotlib.lines as mlines
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter
from PIL import Image


SCRIPT_DIR = Path(__file__).resolve().parent
STORY_DIR = SCRIPT_DIR.parent
PROJECT_ROOT = STORY_DIR.parents[3]
TABLE_DIR = STORY_DIR / "05_tables"
FIGURE_DIR = STORY_DIR / "03_main_figure_candidates"

FIGURE_BASENAME = "Figure_2_candidate_pool_CD34_HSPC_context_typography_20260726_v7"
PAYLOAD_PATH = TABLE_DIR / f"{FIGURE_BASENAME}_source_payload.json"
WORKBOOK_EXTRACT_PATH = TABLE_DIR / "Figure_2_candidate_pool_CD34_HSPC_context_20260716_v1_frozen_workbook_extract.json"

SUPPLEMENTARY_BOOK = TABLE_DIR / "Supplementary_Tables_S1-S6_submission_candidate_20260716_v3b.xlsx"

CD34_RESULT_DIR = (
    PROJECT_ROOT
    / "02_项目总归档_拍板版"
    / "04_前置方法与算法"
    / "方案前置方法"
    / "02_GEO_WGCNA_原始数据与设计"
    / "04_GSE247531_CD34_HSPC分析结果"
)
CD34_EXPRESSION = CD34_RESULT_DIR / "04_intersection126_CD34_log1pCP10K_by_sample.tsv"
CD34_DETECTION = CD34_RESULT_DIR / "05_intersection126_CD34_detection_fraction_by_sample.tsv"

MARKER_DIR = PROJECT_ROOT / "手动下载图片"
MARKER_EXPRESSION = MARKER_DIR / "CD34_marker_log1pCP10K_by_sample.tsv"
MARKER_GROUP_SUMMARY = MARKER_DIR / "CD34_marker_group_summary.tsv"
CD34_SAMPLE_QC = MARKER_DIR / "CD34_sample_QC_by_suffix.tsv"

BULK_RESULT_DIR = (
    PROJECT_ROOT
    / "02_项目总归档_拍板版"
    / "04_前置方法与算法"
    / "方案前置方法"
    / "02_GEO_WGCNA_原始数据与设计"
    / "03_GSE165870_bulk分析结果"
)
BULK_CONTEXT = BULK_RESULT_DIR / "13_intersection126_mapped_to_GSE165870_DESeq2.csv"
BULK_METADATA = BULK_RESULT_DIR / "02_GSE165870_sample_metadata.csv"


BLUE = "#315B8A"
BLUE_LIGHT = "#AFC8E4"
TEAL = "#3A8F8B"
TEAL_LIGHT = "#C8E1DD"
GOLD = "#D7A844"
GOLD_LIGHT = "#F3E3B8"
ORANGE = "#D9822B"
GREY = "#A7ADB4"
GREY_LIGHT = "#E8EAED"
GREY_DARK = "#4A4F55"
BLACK = "#202428"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def rel(path: Path) -> str:
    return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()


def require_files(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Required Figure 2 source files are missing: " + " | ".join(missing))


def records(frame: pd.DataFrame) -> list[dict]:
    clean = frame.astype(object).where(pd.notna(frame), None)
    return clean.to_dict(orient="records")


def expand_suffix_map(metadata: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for row in metadata.to_dict(orient="records"):
        suffixes = [item.strip() for item in str(row["input_count_suffixes"]).split(";") if item.strip()]
        for suffix in suffixes:
            rows.append(
                {
                    "count_suffix": suffix,
                    "subject_timepoint_id": row["subject_timepoint_id"],
                    "subject": row["subject"],
                    "disease": row["disease"],
                    "timepoint": row["timepoint"],
                    "group": row["group"],
                }
            )
    suffix_map = pd.DataFrame(rows)
    if suffix_map["count_suffix"].duplicated().any():
        raise AssertionError("A CD34 input suffix maps to more than one collapsed profile.")
    if len(suffix_map) != 48 or suffix_map["subject_timepoint_id"].nunique() != 42:
        raise AssertionError("Expected 48 input profiles collapsed to 42 subject-by-timepoint profiles.")
    return suffix_map


def reconstruct_collapsed_candidate_data(
    strict_metadata: pd.DataFrame,
    suffix_map: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(CD34_EXPRESSION, sep="\t")
    det = pd.read_csv(CD34_DETECTION, sep="\t")

    expected_pairs = 48 * 126
    if len(expr) != expected_pairs or len(det) != expected_pairs:
        raise AssertionError("Unexpected CD34 candidate-expression source dimensions.")

    expr["reconstructed_count"] = np.rint(
        np.expm1(expr["log1p_cp10k"].astype(float)) * expr["total_umi"].astype(float) / 10000.0
    ).astype(np.int64)
    det["reconstructed_detected_cells"] = np.rint(
        det["detection_fraction"].astype(float) * det["n_cells"].astype(float)
    ).astype(np.int64)

    expr = expr.merge(suffix_map, on="count_suffix", how="left", validate="many_to_one", suffixes=("", "_strict"))
    det = det.merge(suffix_map, on="count_suffix", how="left", validate="many_to_one", suffixes=("", "_strict"))
    if expr["subject_timepoint_id"].isna().any() or det["subject_timepoint_id"].isna().any():
        raise AssertionError("A CD34 input profile is missing from the strict 48-to-42 collapse map.")

    totals = strict_metadata[
        ["subject_timepoint_id", "subject", "disease", "timepoint", "group", "n_cells_sum", "total_umi_from_matrix_sum"]
    ].copy()

    expr42 = (
        expr.groupby(["subject_timepoint_id", "gene"], as_index=False)["reconstructed_count"]
        .sum()
        .merge(totals, on="subject_timepoint_id", how="left", validate="many_to_one")
    )
    expr42["log1p_cp10k"] = np.log1p(
        expr42["reconstructed_count"] / expr42["total_umi_from_matrix_sum"] * 10000.0
    )

    det42 = (
        det.groupby(["subject_timepoint_id", "gene"], as_index=False)["reconstructed_detected_cells"]
        .sum()
        .merge(totals[["subject_timepoint_id", "n_cells_sum"]], on="subject_timepoint_id", how="left", validate="many_to_one")
    )
    det42["detection_fraction"] = det42["reconstructed_detected_cells"] / det42["n_cells_sum"]

    collapsed = expr42.merge(
        det42[["subject_timepoint_id", "gene", "reconstructed_detected_cells", "detection_fraction"]],
        on=["subject_timepoint_id", "gene"],
        how="left",
        validate="one_to_one",
    )
    if len(collapsed) != 42 * 126:
        raise AssertionError("Collapsed CD34 candidate matrix is not 42 profiles by 126 genes.")

    summary = (
        collapsed.groupby("gene", as_index=False)
        .agg(
            n_subject_timepoint_profiles=("subject_timepoint_id", "nunique"),
            mean_log1p_cp10k=("log1p_cp10k", "mean"),
            median_log1p_cp10k=("log1p_cp10k", "median"),
            mean_detection_fraction=("detection_fraction", "mean"),
            median_detection_fraction=("detection_fraction", "median"),
            maximum_detection_fraction=("detection_fraction", "max"),
        )
        .sort_values(["mean_detection_fraction", "gene"], ascending=[False, True])
        .reset_index(drop=True)
    )
    summary["detection_rank_descending"] = np.arange(1, len(summary) + 1)
    return collapsed, summary


def reconstruct_marker_correlations(
    candidate42: pd.DataFrame,
    strict_metadata: pd.DataFrame,
    suffix_map: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    marker = pd.read_csv(MARKER_EXPRESSION, sep="\t")
    marker_map_source = pd.read_csv(MARKER_GROUP_SUMMARY, sep="\t")
    marker_map = marker_map_source[["gene", "panel"]].drop_duplicates()
    sample_qc = pd.read_csv(CD34_SAMPLE_QC, sep="\t")[["count_suffix", "total_umi"]]

    marker = marker.merge(sample_qc, on="count_suffix", how="left", validate="many_to_one")
    marker["reconstructed_count"] = np.rint(
        np.expm1(marker["log1p_cp10k"].astype(float)) * marker["total_umi"].astype(float) / 10000.0
    ).astype(np.int64)
    marker = marker.merge(suffix_map, on="count_suffix", how="left", validate="many_to_one", suffixes=("", "_strict"))
    if marker["subject_timepoint_id"].isna().any():
        raise AssertionError("A marker-panel input profile is missing from the strict collapse map.")

    totals = strict_metadata[["subject_timepoint_id", "total_umi_from_matrix_sum"]]
    marker42 = (
        marker.groupby(["subject_timepoint_id", "gene"], as_index=False)["reconstructed_count"]
        .sum()
        .merge(totals, on="subject_timepoint_id", how="left", validate="many_to_one")
        .merge(marker_map, on="gene", how="left", validate="many_to_one")
    )
    marker42["log1p_cp10k"] = np.log1p(
        marker42["reconstructed_count"] / marker42["total_umi_from_matrix_sum"] * 10000.0
    )
    if marker42["panel"].isna().any():
        raise AssertionError("A marker gene lacks its predefined program annotation.")

    panel42 = (
        marker42.groupby(["subject_timepoint_id", "panel"], as_index=False)
        .agg(panel_score=("log1p_cp10k", "mean"), n_marker_genes=("gene", "nunique"))
    )

    candidate_wide = candidate42.pivot(index="subject_timepoint_id", columns="gene", values="log1p_cp10k")
    panel_wide = panel42.pivot(index="subject_timepoint_id", columns="panel", values="panel_score")
    common_profiles = candidate_wide.index.intersection(panel_wide.index)
    if len(common_profiles) != 42:
        raise AssertionError("Marker-correlation input is not the frozen set of 42 profiles.")

    correlation_rows: list[dict] = []
    for panel in sorted(panel_wide.columns):
        y = panel_wide.loc[common_profiles, panel]
        for gene in sorted(candidate_wide.columns):
            x = candidate_wide.loc[common_profiles, gene]
            if x.nunique(dropna=True) < 2 or y.nunique(dropna=True) < 2:
                rho = np.nan
            else:
                rho = x.rank(method="average").corr(y.rank(method="average"), method="pearson")
            correlation_rows.append(
                {
                    "gene": gene,
                    "panel": panel,
                    "spearman_rho": float(rho),
                    "n_subject_timepoint_profiles": int(len(common_profiles)),
                }
            )
    correlations = pd.DataFrame(correlation_rows)
    correlations["interpretation_boundary"] = (
        "Descriptive profile-level Spearman association after technical-repeat collapse; repeated subjects are not modeled and no inferential P value is used."
    )
    return marker42, panel42, correlations


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.04, label, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="left", va="bottom")


def rounded_box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    facecolor: str,
    edgecolor: str,
    title: str,
    number: str,
    note: str | None = None,
    number_color: str = BLACK,
) -> None:
    x, y = xy
    box = patches.FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.025",
        linewidth=0.8,
        edgecolor=edgecolor,
        facecolor=facecolor,
    )
    ax.add_patch(box)
    ax.text(x + width / 2, y + height * 0.73, title, ha="center", va="center", fontsize=6.0, color=GREY_DARK, linespacing=0.95)
    ax.text(x + width / 2, y + height * 0.39, number, ha="center", va="center", fontsize=11.5, fontweight="bold", color=number_color)
    if note:
        ax.text(x + width / 2, y + height * 0.12, note, ha="center", va="center", fontsize=4.9, color=GREY_DARK)


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = GREY_DARK, dashed: bool = False) -> None:
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops={
            "arrowstyle": "-|>",
            "lw": 0.8,
            "color": color,
            "linestyle": "--" if dashed else "-",
            "mutation_scale": 8,
        },
    )


def plot_panel_a(ax: plt.Axes) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    add_panel_label(ax, "A")
    ax.set_title("Entity parsing and exact candidate-pool construction", loc="left", fontsize=7.4, fontweight="bold", pad=5)

    rounded_box(ax, (0.04, 0.72), 0.34, 0.20, GREY_LIGHT, GREY, "Compound-associated\ntarget entities", "481")
    rounded_box(ax, (0.02, 0.41), 0.32, 0.21, BLUE_LIGHT, BLUE, "Single-gene symbol\nentities", "472", "346 drug-only + 126 shared", BLUE)
    rounded_box(ax, (0.39, 0.41), 0.25, 0.21, GREY_LIGHT, GREY, "Composite target\nentities", "9", "retained outside exact matching", GREY_DARK)
    rounded_box(ax, (0.68, 0.72), 0.30, 0.20, TEAL_LIGHT, TEAL, "AA resource\ngenes", "1,529", "1,403 AA-only + 126 shared", TEAL)
    rounded_box(ax, (0.42, 0.10), 0.41, 0.23, GOLD_LIGHT, GOLD, "Exact symbol\nintersection", "126", "candidate target pool", GOLD)

    ax.plot([0.21, 0.21], [0.72, 0.67], color=GREY_DARK, lw=0.8)
    arrow(ax, (0.21, 0.67), (0.18, 0.62))
    arrow(ax, (0.21, 0.67), (0.50, 0.62))
    arrow(ax, (0.18, 0.41), (0.50, 0.32), BLUE)
    arrow(ax, (0.83, 0.72), (0.73, 0.33), TEAL)
    ax.text(0.50, 0.015, "Only the 472 single-gene entities enter the set operation; box area is not quantitative.", ha="center", fontsize=4.8, color=GREY_DARK)


def plot_panel_b(ax: plt.Axes, candidate_summary: pd.DataFrame) -> None:
    add_panel_label(ax, "B")
    ranked = candidate_summary.sort_values(["mean_detection_fraction", "gene"]).reset_index(drop=True)
    ranked["rank_ascending"] = np.arange(1, len(ranked) + 1)
    ax.plot(ranked["rank_ascending"], ranked["mean_detection_fraction"], color=BLUE_LIGHT, lw=1.0, zorder=1)
    ax.scatter(
        ranked["rank_ascending"],
        ranked["mean_detection_fraction"],
        s=10,
        color=BLUE,
        alpha=0.72,
        linewidth=0,
        zorder=2,
    )
    median = candidate_summary["mean_detection_fraction"].median()
    ax.axhline(median, color=GOLD, lw=1.0, ls="--")
    ax.text(4, median + 0.025, f"median {median:.1%}", fontsize=5.7, color=GOLD, va="bottom")

    top = ranked.nlargest(4, "mean_detection_fraction")
    offsets = {"MIF": (-2, 0.034), "HSP90AB1": (-22, 0.01), "HSP90AA1": (-10, -0.055), "DUT": (-5, 0.035)}
    for row in top.itertuples(index=False):
        dx, dy = offsets.get(row.gene, (-5, 0.025))
        ax.annotate(
            rf"$\mathit{{{row.gene}}}$",
            xy=(row.rank_ascending, row.mean_detection_fraction),
            xytext=(row.rank_ascending + dx, row.mean_detection_fraction + dy),
            fontsize=5.2,
            color=GREY_DARK,
            arrowprops={"arrowstyle": "-", "lw": 0.45, "color": GREY},
        )

    detected = int((candidate_summary["mean_detection_fraction"] > 0).sum())
    ax.text(
        0.03,
        0.94,
        f"126/126 mapped\n{detected}/126 detected\n42 profiles; 23 subjects",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.7,
        color=GREY_DARK,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": GREY_LIGHT, "linewidth": 0.6},
    )
    ax.set_title(r"CD34$^{+}$ HSPC detection coverage", loc="left", fontsize=7.4, fontweight="bold", pad=5)
    ax.set_xlabel("Candidate genes ranked by mean detected-cell fraction", fontsize=6.2)
    ax.set_ylabel("Mean detected-cell fraction", fontsize=6.2)
    ax.set_xlim(0, 128)
    ax.set_ylim(-0.02, max(0.9, ranked["mean_detection_fraction"].max() + 0.08))
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0, decimals=0))
    ax.set_xticks([1, 25, 50, 75, 100, 126])
    ax.tick_params(labelsize=5.5, length=2)
    ax.spines[["top", "right"]].set_visible(False)


PANEL_LABELS = {
    "HSPC_identity": "HSPC identity",
    "HSPC_injury_stress": "HSPC injury / stress",
    "interferon_inflammation": "Interferon / inflammation",
    "hematopoietic_support_context": "Hematopoietic support",
    "liver_THPO_MPL_response": "THPO–MPL response",
}


def plot_panel_c(ax: plt.Axes, correlations: pd.DataFrame) -> None:
    add_panel_label(ax, "C")
    order = [
        "HSPC_identity",
        "HSPC_injury_stress",
        "interferon_inflammation",
        "hematopoietic_support_context",
        "liver_THPO_MPL_response",
    ]
    y_positions = np.arange(len(order))[::-1]
    rng = np.random.default_rng(20260716)
    for panel, y in zip(order, y_positions):
        subset = correlations.loc[correlations["panel"] == panel].copy()
        jitter = rng.uniform(-0.15, 0.15, size=len(subset))
        ax.scatter(
            subset["spearman_rho"],
            np.full(len(subset), y) + jitter,
            s=8,
            color=BLUE,
            alpha=0.28,
            linewidth=0,
            rasterized=False,
        )
        median = subset["spearman_rho"].median()
        ax.scatter([median], [y], marker="D", s=22, color=GOLD, edgecolor="white", linewidth=0.5, zorder=4)
        strongest = subset.nlargest(1, "spearman_rho").iloc[0]
        ax.scatter([strongest["spearman_rho"]], [y], s=20, color=TEAL, edgecolor="white", linewidth=0.5, zorder=5)
        ax.text(
            min(strongest["spearman_rho"] + 0.035, 0.87),
            y + 0.19,
            rf"$\mathit{{{strongest['gene']}}}$  {strongest['spearman_rho']:.2f}",
            fontsize=5.2,
            color=TEAL,
            ha="left",
            va="center",
        )

    ax.axvline(0, color=GREY, lw=0.7)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([PANEL_LABELS[p] for p in order], fontsize=5.6)
    ax.set_xlim(-1, 1)
    ax.set_ylim(-0.55, len(order) - 0.35)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_xlabel("Spearman ρ across 42 subject×timepoint profiles", fontsize=6.2)
    ax.set_title("Association with predefined HSPC marker programs", loc="left", fontsize=7.4, fontweight="bold", pad=5)
    ax.tick_params(axis="x", labelsize=5.5, length=2)
    ax.tick_params(axis="y", length=0)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.text(0.01, 0.03, "dots: candidates   diamonds: median ρ   teal labels: highest positive ρ", transform=ax.transAxes, fontsize=4.8, color=GREY_DARK, ha="left")


def plot_panel_d(ax: plt.Axes, bulk: pd.DataFrame) -> None:
    add_panel_label(ax, "D")
    plotted = bulk.loc[bulk["baseMean"].notna() & bulk["log2FoldChange"].notna()].copy()
    plotted["log10_baseMean_plus1"] = np.log10(plotted["baseMean"].astype(float) + 1.0)
    plotted["nominal_p_lt_0_05"] = plotted["pvalue"].astype(float) < 0.05
    plotted["p_missing"] = plotted["pvalue"].isna()

    other = plotted.loc[~plotted["nominal_p_lt_0_05"] & ~plotted["p_missing"]]
    nominal = plotted.loc[plotted["nominal_p_lt_0_05"]]
    missing = plotted.loc[plotted["p_missing"]]
    ax.scatter(other["log10_baseMean_plus1"], other["log2FoldChange"], s=12, color=GREY, alpha=0.65, linewidth=0)
    ax.scatter(
        missing["log10_baseMean_plus1"],
        missing["log2FoldChange"],
        s=14,
        facecolor="white",
        edgecolor=GREY,
        linewidth=0.6,
    )
    ax.scatter(
        nominal["log10_baseMean_plus1"],
        nominal["log2FoldChange"],
        s=18,
        color=ORANGE,
        edgecolor="white",
        linewidth=0.4,
        zorder=3,
    )
    ax.axhline(0, color=GREY_DARK, lw=0.7)
    ax.axhline(1, color=GREY, lw=0.7, ls="--")
    ax.axhline(-1, color=GREY, lw=0.7, ls="--")

    label_rows = nominal.nsmallest(5, "pvalue")
    annotation_offsets = {
        "IGF1R": (-0.18, -0.55),
        "CA2": (-0.26, 0.55),
        "PIK3CB": (-0.65, -0.15),
        "TEK": (-0.18, 0.55),
        "BCL2L1": (0.12, 0.15),
    }
    for row in label_rows.itertuples(index=False):
        dx, dy = annotation_offsets.get(row.GeneSymbol, (0.05, 0.35))
        ax.annotate(
            rf"$\mathit{{{row.GeneSymbol}}}$",
            xy=(row.log10_baseMean_plus1, row.log2FoldChange),
            xytext=(row.log10_baseMean_plus1 + dx, row.log2FoldChange + dy),
            fontsize=5.0,
            color=GREY_DARK,
            arrowprops={"arrowstyle": "-", "lw": 0.4, "color": GREY},
        )

    fdr_count = int((bulk["padj"].astype(float) < 0.05).sum())
    expressed_count = int(bulk["expressed_in_GSE165870"].astype(bool).sum())
    ax.text(
        0.03,
        0.96,
        f"Healthy n=3; AA n=6\nExpressed: {expressed_count}/126\nFDR < 0.05: {fdr_count}/126",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.6,
        color=GREY_DARK,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": GREY_LIGHT, "linewidth": 0.6},
    )
    ax.set_title(r"Independent Lin$^{-}$CD34$^{+}$ bulk context", loc="left", fontsize=7.4, fontweight="bold", pad=5)
    ax.set_xlabel("log10(baseMean + 1)", fontsize=6.2)
    ax.set_ylabel("log2 fold change (AA / healthy)", fontsize=6.2)
    ax.tick_params(labelsize=5.5, length=2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0.15, 3.85)
    ax.set_ylim(-6.0, 4.6)
    legend = [
        mlines.Line2D([], [], marker="o", ls="", color=ORANGE, markersize=4, label=r"Nominal $P$ < 0.05"),
        mlines.Line2D([], [], marker="o", ls="", color=GREY, markersize=4, label="Other tested"),
        mlines.Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor=GREY, markersize=4, label=r"$P$ unavailable"),
    ]
    ax.legend(handles=legend, loc="lower right", fontsize=4.8, frameon=False, ncol=1, handletextpad=0.35, borderaxespad=0.3)


def main() -> None:
    source_files = [
        SUPPLEMENTARY_BOOK,
        WORKBOOK_EXTRACT_PATH,
        CD34_EXPRESSION,
        CD34_DETECTION,
        MARKER_EXPRESSION,
        MARKER_GROUP_SUMMARY,
        CD34_SAMPLE_QC,
        BULK_CONTEXT,
        BULK_METADATA,
    ]
    require_files(source_files)

    frozen = json.loads(WORKBOOK_EXTRACT_PATH.read_text(encoding="utf-8"))
    entity_records = pd.DataFrame(frozen["S4_drug_targets_481"])
    aa_genes = pd.DataFrame(frozen["S4_AA_targets_1529"])
    candidate_pool = pd.DataFrame(frozen["S4_candidate_pool_126"])
    strict_metadata = pd.DataFrame(frozen["S6_92_metadata"])
    priority = pd.DataFrame(frozen["S5_priority_candidates"])
    for column in [
        "n_input_profiles",
        "n_technical_repeat_flagged_profiles",
        "n_cells_sum",
        "total_umi_metadata_sum",
        "total_umi_from_matrix_sum",
        "total_umi_delta_sum",
    ]:
        strict_metadata[column] = pd.to_numeric(strict_metadata[column], errors="raise")
    priority["rank_no_docking_formal"] = pd.to_numeric(priority["rank_no_docking_formal"], errors="raise")

    single_entities = entity_records.loc[entity_records["entity_type"] == "single_gene_symbol"].copy()
    composite_entities = entity_records.loc[entity_records["entity_type"] == "multi_gene_target_entity"].copy()
    drug_gene_set = set(single_entities["source_target_entity"].astype(str))
    aa_gene_set = set(aa_genes["GeneSymbol"].astype(str))
    candidate_set = set(candidate_pool["GeneSymbol"].astype(str))
    exact_intersection = drug_gene_set & aa_gene_set

    if (len(entity_records), len(single_entities), len(composite_entities)) != (481, 472, 9):
        raise AssertionError("The frozen 481 = 472 + 9 target-entity model did not validate.")
    if len(aa_gene_set) != 1529 or len(candidate_set) != 126 or exact_intersection != candidate_set:
        raise AssertionError("The frozen 472-by-1,529 exact intersection did not reproduce 126 candidates.")
    if composite_entities["candidate_pool_exact_match"].astype(bool).any():
        raise AssertionError("A retained composite entity was incorrectly counted in the candidate pool.")

    top30_annotation = priority[["GeneSymbol", "rank_no_docking_formal", "core10_role"]].copy()
    candidate_annotated = candidate_pool.merge(top30_annotation, on="GeneSymbol", how="left", validate="one_to_one")
    candidate_annotated["candidate_pool_role"] = "exact_single_gene_symbol_intersection"

    suffix_map = expand_suffix_map(strict_metadata)
    candidate42, candidate_summary = reconstruct_collapsed_candidate_data(strict_metadata, suffix_map)
    if set(candidate_summary["gene"]) != candidate_set:
        raise AssertionError("The CD34 expression summary is not an exact mapping of the 126-gene pool.")

    marker42, panel42, correlations = reconstruct_marker_correlations(candidate42, strict_metadata, suffix_map)

    bulk = pd.read_csv(BULK_CONTEXT)
    bulk_metadata = pd.read_csv(BULK_METADATA)
    bulk_metadata = bulk_metadata.drop(columns=[column for column in bulk_metadata.columns if str(column).startswith("Unnamed:")])
    if set(bulk["GeneSymbol"].astype(str)) != candidate_set:
        raise AssertionError("The GSE165870 candidate context is not an exact 126-gene mapping table.")
    if bulk_metadata["condition"].value_counts().to_dict() != {"AA": 6, "Healthy": 3}:
        raise AssertionError("Unexpected GSE165870 group counts.")

    entity_summary = pd.DataFrame(
        [
            {"metric": "compound_associated_target_entities", "count": 481, "definition": "Deduplicated source target entities"},
            {"metric": "single_gene_symbol_entities", "count": 472, "definition": "Eligible for exact gene-symbol matching"},
            {"metric": "retained_composite_target_entities", "count": 9, "definition": "Retained without decomposition; excluded from exact intersection"},
            {"metric": "AA_resource_genes", "count": 1529, "definition": "Deduplicated AA resource gene symbols"},
            {"metric": "candidate_target_pool", "count": 126, "definition": "Exact intersection of 472 single-gene entities and 1,529 AA resource genes"},
            {"metric": "drug_only_single_gene_entities", "count": 346, "definition": "472 minus 126"},
            {"metric": "AA_only_resource_genes", "count": 1403, "definition": "1,529 minus 126"},
        ]
    )

    figure_summary = pd.DataFrame(
        [
            {"metric": "candidate_genes_mapped_to_CD34", "value": 126},
            {"metric": "candidate_genes_detected_in_CD34", "value": int((candidate_summary["mean_detection_fraction"] > 0).sum())},
            {"metric": "CD34_subject_timepoint_profiles", "value": 42},
            {"metric": "CD34_subjects", "value": int(strict_metadata["subject"].nunique())},
            {"metric": "CD34_cells", "value": int(strict_metadata["n_cells_sum"].sum())},
            {"metric": "GSE165870_healthy_samples", "value": 3},
            {"metric": "GSE165870_AA_samples", "value": 6},
            {"metric": "GSE165870_candidates_expressed", "value": int(bulk["expressed_in_GSE165870"].astype(bool).sum())},
            {"metric": "GSE165870_candidate_FDR_lt_0_05", "value": int((bulk["padj"].astype(float) < 0.05).sum())},
        ]
    )

    source_index = pd.DataFrame(
        [
            {
                "source_role": role,
                "project_relative_path": rel(path),
                "sha256": sha256(path),
            }
            for role, path in [
                ("Frozen supplementary source workbook", SUPPLEMENTARY_BOOK),
                ("Candidate expression by input profile", CD34_EXPRESSION),
                ("Candidate detection by input profile", CD34_DETECTION),
                ("Marker-gene expression by input profile", MARKER_EXPRESSION),
                ("Marker program dictionary", MARKER_GROUP_SUMMARY),
                ("CD34 input-profile QC", CD34_SAMPLE_QC),
                ("Independent bulk candidate context", BULK_CONTEXT),
                ("Independent bulk sample metadata", BULK_METADATA),
            ]
        ]
    )

    payload = {
        "figure": FIGURE_BASENAME,
        "figure_contract": {
            "core_conclusion": (
                "An exact 126-gene candidate space arises only from standardized single-gene entities and is subsequently contextualized, not redefined, by CD34⁺ HSPC expression and independent bulk data."
            ),
            "archetype": "asymmetric mixed-modality figure",
            "backend": "Python",
            "final_size_mm": {"width": 183, "height": 150},
            "statistics_boundary": (
                "CD34 marker correlations are descriptive after technical-repeat collapse and carry no inferential P value; GSE165870 differential expression uses DESeq2 and BH adjustment."
            ),
        },
        "entity_summary": records(entity_summary),
        "entity_records": records(entity_records),
        "AA_resource_genes": records(aa_genes),
        "candidate_pool_126": records(candidate_annotated),
        "CD34_strict_metadata_42": records(strict_metadata),
        "CD34_candidate_expression_42": records(candidate42),
        "CD34_candidate_summary_126": records(candidate_summary),
        "CD34_marker_gene_expression_42": records(marker42),
        "CD34_marker_program_scores_42": records(panel42),
        "CD34_marker_correlations_126x5": records(correlations),
        "GSE165870_metadata": records(bulk_metadata),
        "GSE165870_candidate_context": records(bulk),
        "figure_summary": records(figure_summary),
        "source_index": records(source_index),
    }
    PAYLOAD_PATH.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    width_in = 183.0 / 25.4
    height_in = 150.0 / 25.4
    fig = plt.figure(figsize=(width_in, height_in), facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.07, 1.0],
        height_ratios=[0.98, 1.02],
        left=0.12,
        right=0.99,
        bottom=0.10,
        top=0.965,
        wspace=0.27,
        hspace=0.38,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    plot_panel_a(ax_a)
    plot_panel_b(ax_b, candidate_summary)
    plot_panel_c(ax_c, correlations)
    plot_panel_d(ax_d, bulk)

    output_base = FIGURE_DIR / FIGURE_BASENAME
    fig.savefig(output_base.with_suffix(".svg"), facecolor="white")
    fig.savefig(output_base.with_suffix(".pdf"), facecolor="white")
    png_path = output_base.with_suffix(".png")
    fig.savefig(png_path, dpi=600, facecolor="white")
    plt.close(fig)
    with Image.open(png_path) as png_image:
        png_image.convert("RGB").save(png_path, dpi=(600, 600), optimize=True)

    print(json.dumps({
        "figure_base": str(output_base),
        "payload": str(PAYLOAD_PATH),
        "candidate_count": len(candidate_set),
        "CD34_profiles": len(strict_metadata),
        "CD34_subjects": int(strict_metadata["subject"].nunique()),
        "CD34_cells": int(strict_metadata["n_cells_sum"].sum()),
        "bulk_FDR_lt_0_05": int((bulk["padj"].astype(float) < 0.05).sum()),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
