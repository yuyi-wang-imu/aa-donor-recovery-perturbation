from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
DATA = Path(
    os.environ.get("AA_GENEFORMER_DATA_DIR", ROOT / "derived_data" / "geneformer")
)
OUT = Path(
    os.environ.get("AA_GENEFORMER_FIGURE_OUT", ROOT / "outputs" / "geneformer_figures")
)
OUT.mkdir(parents=True, exist_ok=True)

AUDIT = DATA / "AA_Geneformer_BIDIRECTIONAL_AUDIT_SINGLE_20260802_v1"
EXT = DATA / "AA_Geneformer_EXTENSION_TABLES_20260802_v1"
POS = DATA / "AA_Geneformer_POSCTRL_MATCHED_20260802_v2"

FILES = {
    "candidate": AUDIT / "candidate_bidirectional_summary.tsv",
    "candidate_boot": AUDIT / "candidate_donor_bootstrap.tsv",
    "donor_bidir": AUDIT / "donor_paired_bidirectional_effects.tsv",
    "baselines": EXT / "Table_C1_simple_baselines.tsv",
    "lodo": EXT / "Table_C2_LODO_by_donor.tsv",
    "ablation": EXT / "Table_D_program_ablation.tsv",
    "matching": DATA / "AA_Geneformer_REMATCH_20260802_v1_diagnostics.tsv",
    "matched_null": DATA / "AA_Geneformer_REMATCH_20260802_v1_candidate_matched_null.tsv",
    "state_boot": DATA / "AA_Geneformer_STATE_SPECIFIC_20260802_v1_state_gene_bootstrap.tsv",
    "positive": POS / "positive_control_gene_summary.tsv",
    "cross_model": DATA / "AA_Geneformer_FULL_20260802_v1_geneformer_sctenifold_consistency.tsv",
    "recovery": DATA / "AA_candidate_set_directional_recovery_extension_MVP_20260802_v1_donor_scores.csv",
    "external_145668": DATA / "AA_GSE145668_external_paired_recovery_MVP_20260802_v1_donor_scores.csv",
    "external_145668_json": DATA / "AA_GSE145668_external_paired_recovery_MVP_20260802_v1.json",
    "external_165870": DATA / "AA_GSE165870_candidate_direction_replication_MVP_20260802_v1_gene_effects.csv",
    "external_165870_json": DATA / "AA_GSE165870_candidate_direction_replication_MVP_20260802_v1.json",
}

for key, path in FILES.items():
    if not path.exists():
        raise FileNotFoundError(f"Missing input for {key}: {path}")


available_fonts = {font.name for font in font_manager.fontManager.ttflist}
font_family = next(
    (name for name in ("Arial", "Liberation Sans", "DejaVu Sans") if name in available_fonts),
    "sans-serif",
)

mpl.rcParams.update(
    {
        "font.family": font_family,
        "font.size": 6.5,
        "axes.titlesize": 7,
        "axes.labelsize": 6.5,
        "xtick.labelsize": 5.5,
        "ytick.labelsize": 5.5,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "legend.fontsize": 5.5,
        "lines.linewidth": 0.8,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.unicode_minus": False,
    }
)

BLUE = "#3B6FB6"
ORANGE = "#D9781F"
TEAL = "#2A9D8F"
PURPLE = "#7A5AA6"
RED = "#B84A4A"
GRAY = "#7A8088"
LIGHT = "#E8ECF1"
DARK = "#20242A"

PANEL_SIZE = (3.504, 2.72)  # 89 mm wide; convenient for later 2-column assembly.
TITLE_Y = 0.955
LETTER_X = 0.035
TITLE_X = 0.115


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() == ".tsv" else ","
    return pd.read_csv(path, sep=sep)


candidate = read_table(FILES["candidate"])
candidate_boot = read_table(FILES["candidate_boot"])
donor_bidir = read_table(FILES["donor_bidir"])
baselines = read_table(FILES["baselines"])
lodo = read_table(FILES["lodo"])
ablation = read_table(FILES["ablation"])
matching = read_table(FILES["matching"])
matched_null = read_table(FILES["matched_null"])
state_boot = read_table(FILES["state_boot"])
positive = read_table(FILES["positive"])
cross_model = read_table(FILES["cross_model"])
recovery = read_table(FILES["recovery"])
external_145668 = read_table(FILES["external_145668"])
external_165870 = read_table(FILES["external_165870"])
with FILES["external_145668_json"].open("r", encoding="utf-8") as handle:
    external_145668_meta = json.load(handle)
with FILES["external_165870_json"].open("r", encoding="utf-8") as handle:
    external_165870_meta = json.load(handle)


def make_panel(letter: str, title: str, left: float = 0.19, bottom: float = 0.18,
               right: float = 0.965, top: float = 0.82):
    fig = plt.figure(figsize=PANEL_SIZE, facecolor="white")
    ax = fig.add_axes([left, bottom, right - left, top - bottom])
    fig.text(LETTER_X, TITLE_Y, letter, ha="left", va="top", fontsize=8,
             fontweight="bold", color=DARK)
    fig.text(TITLE_X, TITLE_Y, title, ha="left", va="top", fontsize=7,
             fontweight="bold", color=DARK)
    return fig, ax


def clean_axes(ax, grid_axis: str | None = "x"):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color="#D9DEE5", linewidth=0.45, zorder=0)
    ax.tick_params(pad=2)


def save_panel(fig, stem: str):
    png = OUT / f"{stem}.png"
    svg = OUT / f"{stem}.svg"
    fig.savefig(png, dpi=600, facecolor="white", bbox_inches=None)
    fig.savefig(svg, facecolor="white", bbox_inches=None)
    plt.close(fig)
    return png, svg


def emphasize_colors(names):
    return [TEAL if x == "TOP2A" else PURPLE if x == "GSK3B" else ORANGE if x == "KIT" else BLUE for x in names]


def panel_9a():
    fig, ax = make_panel("A", "Donor-aware perturbation design", left=0.04, bottom=0.12, right=0.97, top=0.82)
    ax.set_axis_off()
    boxes = [
        (0.03, 0.56, 0.18, 0.24, "Paired donor\ntranscriptomes\n17 donors", BLUE),
        (0.28, 0.56, 0.18, 0.24, "Observed\nrecovery axis\nLODO", TEAL),
        (0.53, 0.56, 0.18, 0.24, "Geneformer\ndelete and\noverexpress", PURPLE),
        (0.78, 0.56, 0.18, 0.24, "Matched\nbackground\ncalibration", ORANGE),
        (0.28, 0.12, 0.18, 0.20, "Donor\nbootstrap", GRAY),
        (0.53, 0.12, 0.18, 0.20, "Cell-state\nspecificity", GRAY),
    ]
    for x, y, w, h, txt, color in boxes:
        patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.018",
                               linewidth=0.8, edgecolor=color, facecolor="white")
        ax.add_patch(patch)
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=6,
                color=DARK, linespacing=1.15)
    for x0, x1 in [(0.21, 0.28), (0.46, 0.53), (0.71, 0.78)]:
        ax.annotate("", xy=(x1, 0.68), xytext=(x0, 0.68),
                    arrowprops=dict(arrowstyle="->", lw=0.8, color=GRAY))
    ax.annotate("", xy=(0.37, 0.32), xytext=(0.37, 0.56),
                arrowprops=dict(arrowstyle="->", lw=0.8, color=GRAY))
    ax.annotate("", xy=(0.62, 0.32), xytext=(0.62, 0.56),
                arrowprops=dict(arrowstyle="->", lw=0.8, color=GRAY))
    ax.text(0.5, 0.01, "Prespecified candidates are evaluated after donor-level aggregation",
            ha="center", va="bottom", fontsize=5.5, color=GRAY)
    return save_panel(fig, "Figure_9A_Donor_Aware_Perturbation_Design_20260803_v1")


def panel_9b():
    df = recovery.loc[recovery["test_id"] == "primary_all10_equal"].drop_duplicates("subject").copy()
    df = df.sort_values("independent_recovery_shift").reset_index(drop=True)
    fig, ax = make_panel("B", "Observed donor-level recovery", left=0.20)
    y = np.arange(len(df))
    vals = df["independent_recovery_shift"].to_numpy()
    ax.hlines(y, 0, vals, color="#BBC6D3", lw=0.7, zorder=1)
    ax.scatter(vals, y, s=20, color=TEAL, edgecolor="white", linewidth=0.35, zorder=3)
    ax.axvline(0, color=DARK, lw=0.7)
    ax.axvline(np.median(vals), color=ORANGE, lw=0.8, ls="--")
    ax.set_yticks(y, df["subject"])
    ax.set_xlabel("Recovery shift toward healthy state")
    ax.set_ylabel("Donor")
    ax.text(0.98, 0.97, f"17/17 positive\nmedian = {np.median(vals):.3f}\nexact $P$ = 1.53 × 10$^{{-5}}$",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.5,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#CCD3DC", linewidth=0.5))
    clean_axes(ax)
    return save_panel(fig, "Figure_9B_Observed_Donor_Level_Recovery_20260803_v1")


def panel_9c():
    df = candidate.copy()
    fig, ax = make_panel("C", "Bidirectional Geneformer perturbation", left=0.18)
    xmin = min(-0.06, df["deletion_recovery_shift"].min() - 0.006)
    xmax = max(0.015, df["deletion_recovery_shift"].max() + 0.006)
    ymin = min(-0.02, df["overexpression_recovery_shift_detected_cells"].min() - 0.006)
    ymax = max(0.045, df["overexpression_recovery_shift_detected_cells"].max() + 0.006)
    ax.axvspan(xmin, 0, ymin=(0 - ymin) / (ymax - ymin), ymax=1, color="#E8F4EF", zorder=0)
    ax.axhline(0, color=GRAY, lw=0.6, ls="--")
    ax.axvline(0, color=GRAY, lw=0.6, ls="--")
    colors = emphasize_colors(df["candidate"])
    ax.scatter(df["deletion_recovery_shift"], df["overexpression_recovery_shift_detected_cells"],
               s=26, c=colors, edgecolor="white", linewidth=0.4, zorder=3)
    for _, row in df.iterrows():
        ax.annotate(row["candidate"],
                    (row["deletion_recovery_shift"], row["overexpression_recovery_shift_detected_cells"]),
                    xytext=(3, 2), textcoords="offset points", fontsize=5.1)
    ax.text(0.03, 0.96, "Expected bidirectional pattern", transform=ax.transAxes,
            ha="left", va="top", fontsize=5.4, color=TEAL)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("Deletion shift toward healthy state")
    ax.set_ylabel("Overexpression shift toward healthy state")
    clean_axes(ax, None)
    return save_panel(fig, "Figure_9C_Bidirectional_Geneformer_Perturbation_20260803_v1")


def panel_9d():
    df = candidate_boot.sort_values("mean_donor_bidirectional_score", ascending=True).reset_index(drop=True)
    fig, ax = make_panel("D", "Donor-level bidirectional effects", left=0.25)
    y = np.arange(len(df))
    x = df["mean_donor_bidirectional_score"].to_numpy()
    lo = df["bootstrap_ci025"].to_numpy()
    hi = df["bootstrap_ci975"].to_numpy()
    colors = emphasize_colors(df["candidate"])
    ax.errorbar(x, y, xerr=np.vstack([x - lo, hi - x]), fmt="none", ecolor="#8993A0",
                elinewidth=0.8, capsize=1.8, zorder=1)
    ax.scatter(x, y, s=24, c=colors, edgecolor="white", linewidth=0.35, zorder=3)
    ax.axvline(0, color=DARK, lw=0.7)
    ax.set_yticks(y, df["candidate"])
    ax.set_xlabel("Mean bidirectional score (95% donor bootstrap CI)")
    clean_axes(ax)
    return save_panel(fig, "Figure_9D_Donor_Level_Bidirectional_Effects_20260803_v1")


def panel_9e():
    df = candidate.sort_values("matched_control_percentile", ascending=True).reset_index(drop=True)
    fig, ax = make_panel("E", "Matched-background calibration", left=0.25)
    y = np.arange(len(df))
    pct = 100 * df["matched_control_percentile"].to_numpy()
    allowed = df["formal_matched_claim_allowed"].astype(bool).to_numpy()
    ax.hlines(y, 0, pct, color="#CBD3DD", lw=0.8)
    ax.scatter(pct[allowed], y[allowed], s=23, color=BLUE, edgecolor="white", linewidth=0.35, zorder=3)
    ax.scatter(pct[~allowed], y[~allowed], s=26, facecolor="white", edgecolor=ORANGE,
               linewidth=0.9, zorder=3)
    ax.axvline(95, color=ORANGE, lw=0.8, ls="--")
    ax.set_xlim(0, 103)
    ax.set_yticks(y, df["candidate"])
    ax.set_xlabel("Percentile within expression-matched background")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE, markeredgecolor="white",
               markersize=4.5, label="Balance criteria met"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=ORANGE,
               markersize=4.5, label="Balance criteria not met"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False, handletextpad=0.4)
    clean_axes(ax)
    return save_panel(fig, "Figure_9E_Matched_Background_Calibration_20260803_v1")


def panel_9f():
    df = state_boot.loc[
        (state_boot["state_class"] == "HSPC-marker-class")
        & (state_boot["metric"] == "bidirectional_recovery_score")
        & (state_boot["gene"].isin(candidate["candidate"]))
    ].copy()
    df = df.sort_values("mean", ascending=True).reset_index(drop=True)
    fig, ax = make_panel("F", "HSPC-state perturbation specificity", left=0.25)
    y = np.arange(len(df))
    x = df["mean"].to_numpy()
    lo = df["ci_low"].to_numpy()
    hi = df["ci_high"].to_numpy()
    ax.errorbar(x, y, xerr=np.vstack([x - lo, hi - x]), fmt="none", ecolor="#8993A0",
                elinewidth=0.8, capsize=1.8)
    ax.scatter(x, y, s=23, c=emphasize_colors(df["gene"]), edgecolor="white", linewidth=0.35, zorder=3)
    ax.axvline(0, color=DARK, lw=0.7)
    ax.set_yticks(y, df["gene"])
    ax.set_xlabel("Mean bidirectional score (95% donor bootstrap CI)")
    ax.text(0.98, 0.03, "HSPC-marker-class cells", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=5.4, color=GRAY)
    clean_axes(ax)
    return save_panel(fig, "Figure_9F_HSPC_State_Perturbation_Specificity_20260803_v1")


def panel_s9a():
    df = baselines.copy()
    labels = {
        "expression_centroid_shift": "Expression centroid",
        "ridge_shift": "Ridge recovery axis",
        "geneformer_observed_recovery_shift": "Geneformer embedding",
    }
    df["label"] = df["baseline"].map(labels)
    fig, ax = make_panel("A", "Recovery-axis baseline comparison", left=0.27)
    y = np.arange(len(df))[::-1]
    vals = df["positive_fraction"].to_numpy()
    ax.barh(y, vals, height=0.48, color=[GRAY, BLUE, PURPLE], edgecolor="none")
    for yi, val, n in zip(y, vals, df["n_donors"]):
        ax.text(min(val + 0.018, 1.04), yi, f"{int(round(val*n))}/{int(n)}", va="center", fontsize=5.5)
    ax.set_yticks(y, df["label"])
    ax.set_xlim(0, 1.08)
    ax.set_xlabel("Fraction of donors shifting toward healthy state")
    ax.text(0.02, 0.03, "Geneformer vs centroid/ridge rank correlations: 0.650 / 0.675",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=5.2, color=GRAY)
    clean_axes(ax)
    return save_panel(fig, "Figure_S9A_Recovery_Axis_Baseline_Comparison_20260803_v1")


def panel_s9b():
    df = matching.sort_values("max_abs_standardized_residual", ascending=True).reset_index(drop=True)
    fig, ax = make_panel("B", "Matched-control covariate balance", left=0.25)
    y = np.arange(len(df))
    vals = df["max_abs_standardized_residual"].to_numpy()
    allowed = df["formal_matched_claim_allowed"].astype(bool).to_numpy()
    ax.barh(y, vals, height=0.48, color=np.where(allowed, BLUE, ORANGE))
    ax.axvline(1.0, color=DARK, lw=0.7, ls="--")
    ax.set_yticks(y, df["candidate"])
    ax.set_xlabel("Maximum absolute standardized residual")
    ax.text(0.98, 0.04, "Reference = 1.0", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.3)
    clean_axes(ax)
    return save_panel(fig, "Figure_S9B_Matched_Control_Covariate_Balance_20260803_v1")


def panel_s9c():
    df = candidate.sort_values("bidirectional_score_min_of_arms", ascending=True).reset_index(drop=True)
    fig, ax = make_panel("C", "Candidate effects versus matched controls", left=0.25)
    y = np.arange(len(df))
    means = df["matched_control_score_mean"].to_numpy()
    sds = df["matched_control_score_sd"].to_numpy()
    scores = df["bidirectional_score_min_of_arms"].to_numpy()
    ax.errorbar(means, y, xerr=sds, fmt="o", ms=3.1, color=GRAY, ecolor="#AAB2BC",
                elinewidth=0.75, capsize=1.6, label="Matched controls")
    ax.scatter(scores, y, s=22, c=emphasize_colors(df["candidate"]), edgecolor="white",
               linewidth=0.35, zorder=3, label="Candidate")
    ax.axvline(0, color=DARK, lw=0.7)
    ax.set_yticks(y, df["candidate"])
    ax.set_xlabel("Bidirectional score")
    ax.legend(loc="lower right", frameon=False)
    clean_axes(ax)
    return save_panel(fig, "Figure_S9C_Candidate_Effects_Versus_Matched_Controls_20260803_v1")


def panel_s9d():
    df = lodo.reset_index(drop=True)
    fig, ax = make_panel("D", "Leave-one-donor-out stability", left=0.16)
    x = np.arange(len(df))
    ax.plot(x, df["spearman_rank_stability"], marker="o", ms=3, color=BLUE, label="Spearman rank")
    ax.plot(x, df["top3_jaccard"], marker="s", ms=2.8, color=ORANGE, label="Top-3 Jaccard")
    ax.set_ylim(0.86, 1.012)
    ax.set_xticks(x, df["left_out_subject"], rotation=55, ha="right")
    ax.set_ylabel("Stability")
    ax.set_xlabel("Left-out donor")
    ax.legend(loc="lower left", frameon=False)
    clean_axes(ax, "y")
    return save_panel(fig, "Figure_S9D_Leave_One_Donor_Out_Stability_20260803_v1")


def panel_s9e():
    pivot = donor_bidir.pivot(index="subject", columns="gene", values="bidirectional_score_min_of_arms")
    order = candidate.sort_values("bidirectional_score_min_of_arms", ascending=False)["candidate"].tolist()
    pivot = pivot.reindex(columns=order).sort_index()
    fig, ax = make_panel("e", "Donor-by-candidate bidirectional effects", left=0.17, bottom=0.27, top=0.82)
    vmax = float(np.nanmax(np.abs(pivot.to_numpy())))
    im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   interpolation="nearest")
    ax.set_xticks(np.arange(len(pivot.columns)), pivot.columns, rotation=55, ha="right")
    ax.set_yticks(np.arange(len(pivot.index)), pivot.index)
    ax.set_xlabel("Candidate gene")
    ax.set_ylabel("Donor")
    # Keep the colour key above the plotting area.  The earlier lower-right
    # placement crossed the rotated CD38/TERT tick labels after the six-panel
    # layout was reduced for the supplementary-information page.
    cax = fig.add_axes([0.75, 0.845, 0.20, 0.025])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.ax.set_title("Bidirectional score", fontsize=5.3, pad=1.5)
    cb.ax.tick_params(labelsize=5, length=2)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    return save_panel(fig, "Figure_S9E_Donor_By_Candidate_Heatmap_20260803_v1")


def panel_s9f():
    df = positive.sort_values("bidirectional_recovery_score", ascending=True).reset_index(drop=True)
    fig, ax = make_panel("F", "Known hematopoietic regulators", left=0.25)
    y = np.arange(len(df))
    good = df["expected_direction_both"].astype(bool).to_numpy()
    vals = df["bidirectional_recovery_score"].to_numpy()
    ax.hlines(y, 0, vals, color="#C8D0DA", lw=0.8)
    ax.scatter(vals[good], y[good], s=24, color=TEAL, edgecolor="white", linewidth=0.35, zorder=3)
    ax.scatter(vals[~good], y[~good], s=24, facecolor="white", edgecolor=GRAY, linewidth=0.8, zorder=3)
    ax.axvline(0, color=DARK, lw=0.7)
    ax.set_yticks(y, df["gene"])
    ax.set_xlabel("Bidirectional recovery score")
    ax.text(0.98, 0.04, "MPL: insufficient baseline coverage", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=5.3, color=GRAY)
    clean_axes(ax)
    return save_panel(fig, "Figure_S9F_Known_Hematopoietic_Regulators_20260803_v1")


def panel_s10a():
    df = ablation.copy()
    labels = {
        "all_programs": "All programs",
        "leave_out_HSPC_identity": "Without HSPC identity",
        "leave_out_cell_cycle_recovery": "Without cell-cycle recovery",
        "leave_out_hematopoietic_support": "Without hematopoietic support",
    }
    df["label"] = df["variant"].map(labels)
    df = df.sort_values("spearman_vs_all_programs", ascending=True).reset_index(drop=True)
    fig, ax = make_panel("A", "Recovery-program ablation", left=0.39)
    y = np.arange(len(df))
    vals = df["spearman_vs_all_programs"].to_numpy()
    colors = [ORANGE if v == "leave_out_cell_cycle_recovery" else BLUE for v in df["variant"]]
    ax.barh(y, vals, height=0.48, color=colors)
    ax.set_yticks(y, df["label"])
    ax.set_xlim(0, 1.04)
    ax.set_xlabel("Spearman correlation with complete recovery axis")
    clean_axes(ax)
    return save_panel(fig, "Figure_S10A_Recovery_Program_Ablation_20260803_v1")


def panel_s10b():
    df = external_145668.copy()
    fig, ax = make_panel("B", "External paired-cohort assessment", left=0.19)
    ax.axhline(0, color="#C9CFD7", lw=0.6)
    ax.axvline(0, color="#C9CFD7", lw=0.6)
    ax.scatter(df["independent_global_recovery_shift"], df["candidate_set_directional_score"],
               s=26, color=BLUE, edgecolor="white", linewidth=0.4)
    for _, row in df.iterrows():
        ax.annotate(row["donor"], (row["independent_global_recovery_shift"], row["candidate_set_directional_score"]),
                    xytext=(3, 2), textcoords="offset points", fontsize=5.2)
    ax.set_xlabel("Observed recovery shift")
    ax.set_ylabel("Candidate-set directional score")
    ax.text(0.98, 0.04,
            f"Spearman ρ = {external_145668_meta['candidate_spearman_rho']:.3f}\npermutation $P$ = {external_145668_meta['exact_permutation_p']:.3f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.4)
    clean_axes(ax, None)
    return save_panel(fig, "Figure_S10B_External_Paired_Cohort_Assessment_20260803_v1")


def panel_s10c():
    df = external_165870.copy()
    fig, ax = make_panel("C", "Cross-cohort disease-direction assessment", left=0.19)
    ax.axhline(0, color="#C9CFD7", lw=0.6)
    ax.axvline(0, color="#C9CFD7", lw=0.6)
    colors = np.where(df["direction_concordant"].astype(bool), TEAL, GRAY)
    ax.scatter(df["GSE247531_healthy_minus_SAA_baseline"], df["GSE165870_healthy_minus_AA_log2FC"],
               s=25, c=colors, edgecolor="white", linewidth=0.35)
    for _, row in df.iterrows():
        ax.annotate(row["gene"],
                    (row["GSE247531_healthy_minus_SAA_baseline"], row["GSE165870_healthy_minus_AA_log2FC"]),
                    xytext=(3, 2), textcoords="offset points", fontsize=5.0)
    ax.set_xlabel("GSE247531: healthy minus SAA baseline")
    ax.set_ylabel("GSE165870: healthy minus AA (log$_2$FC)")
    ax.text(0.98, 0.04,
            f"7/10 directionally concordant\nSpearman ρ = {external_165870_meta['candidate_spearman_rho']:.3f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.4)
    clean_axes(ax, None)
    return save_panel(fig, "Figure_S10C_Cross_Cohort_Disease_Direction_20260803_v1")


def panel_s10d():
    df = cross_model.copy()
    fig, ax = make_panel("D", "Complementary perturbation outputs", left=0.19)
    ax.scatter(df["geneformer_signed_deletion_shift"], df["sctenifold_significant_response_genes"],
               s=25, c=emphasize_colors(df["candidate"]), edgecolor="white", linewidth=0.35)
    for _, row in df.iterrows():
        ax.annotate(row["candidate"],
                    (row["geneformer_signed_deletion_shift"], row["sctenifold_significant_response_genes"]),
                    xytext=(3, 2), textcoords="offset points", fontsize=5.0)
    ax.axvline(0, color="#C9CFD7", lw=0.6)
    ax.set_xlabel("Geneformer deletion shift")
    ax.set_ylabel("scTenifoldKnk response-gene count")
    ax.text(0.98, 0.96, "Spearman ρ = -0.130", transform=ax.transAxes,
            ha="right", va="top", fontsize=5.4)
    clean_axes(ax, None)
    return save_panel(fig, "Figure_S10D_Complementary_Perturbation_Outputs_20260803_v1")


def make_layout_proof(stems, output_stem, ncols=2, bg="#F3F5F7"):
    imgs = [Image.open(OUT / f"{stem}.png").convert("RGB") for stem in stems]
    target_w = 1500
    resized = []
    for img in imgs:
        h = int(round(img.height * target_w / img.width))
        resized.append(img.resize((target_w, h), Image.Resampling.LANCZOS))
    gutter = 70
    outer = 90
    nrows = int(np.ceil(len(resized) / ncols))
    cell_h = max(img.height for img in resized)
    canvas_w = outer * 2 + ncols * target_w + (ncols - 1) * gutter
    canvas_h = outer * 2 + nrows * cell_h + (nrows - 1) * gutter
    canvas = Image.new("RGB", (canvas_w, canvas_h), bg)
    for i, img in enumerate(resized):
        row, col = divmod(i, ncols)
        x = outer + col * (target_w + gutter)
        y = outer + row * (cell_h + gutter)
        canvas.paste(img, (x, y))
    canvas.save(OUT / f"{output_stem}.png", dpi=(300, 300))


def write_source_map(outputs):
    rows = []
    for panel, input_keys in outputs.items():
        for key in input_keys:
            rows.append({"panel": panel, "source_role": key, "source_path": str(FILES[key])})
    pd.DataFrame(rows).to_csv(OUT / "New_Figure_Panel_Source_Map_20260803_v1.tsv", sep="\t", index=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    panel_9a(); panel_9b(); panel_9c(); panel_9d(); panel_9e(); panel_9f()
    panel_s9a(); panel_s9b(); panel_s9c(); panel_s9d(); panel_s9e(); panel_s9f()
    panel_s10a(); panel_s10b(); panel_s10c(); panel_s10d()

    figure9 = [
        "Figure_9A_Donor_Aware_Perturbation_Design_20260803_v1",
        "Figure_9B_Observed_Donor_Level_Recovery_20260803_v1",
        "Figure_9C_Bidirectional_Geneformer_Perturbation_20260803_v1",
        "Figure_9D_Donor_Level_Bidirectional_Effects_20260803_v1",
        "Figure_9E_Matched_Background_Calibration_20260803_v1",
        "Figure_9F_HSPC_State_Perturbation_Specificity_20260803_v1",
    ]
    figures9 = [
        "Figure_S9A_Recovery_Axis_Baseline_Comparison_20260803_v1",
        "Figure_S9B_Matched_Control_Covariate_Balance_20260803_v1",
        "Figure_S9C_Candidate_Effects_Versus_Matched_Controls_20260803_v1",
        "Figure_S9D_Leave_One_Donor_Out_Stability_20260803_v1",
        "Figure_S9E_Donor_By_Candidate_Heatmap_20260803_v1",
        "Figure_S9F_Known_Hematopoietic_Regulators_20260803_v1",
    ]
    figures10 = [
        "Figure_S10A_Recovery_Program_Ablation_20260803_v1",
        "Figure_S10B_External_Paired_Cohort_Assessment_20260803_v1",
        "Figure_S10C_Cross_Cohort_Disease_Direction_20260803_v1",
        "Figure_S10D_Complementary_Perturbation_Outputs_20260803_v1",
    ]
    make_layout_proof(figure9, "Figure_9_Layout_Proof_20260803_v1")
    make_layout_proof(figures9, "Figure_S9_Layout_Proof_20260803_v1")
    make_layout_proof(figures10, "Figure_S10_Layout_Proof_20260803_v1")

    write_source_map(
        {
            "Figure 9A": ["recovery", "candidate"],
            "Figure 9B": ["recovery"],
            "Figure 9C": ["candidate"],
            "Figure 9D": ["candidate_boot"],
            "Figure 9E": ["candidate", "matching"],
            "Figure 9F": ["state_boot"],
            "Figure S9A": ["baselines"],
            "Figure S9B": ["matching"],
            "Figure S9C": ["candidate"],
            "Figure S9D": ["lodo"],
            "Figure S9E": ["donor_bidir"],
            "Figure S9F": ["positive"],
            "Figure S10A": ["ablation"],
            "Figure S10B": ["external_145668", "external_145668_json"],
            "Figure S10C": ["external_165870", "external_165870_json"],
            "Figure S10D": ["cross_model"],
        }
    )
    print(f"Generated {len(figure9) + len(figures9) + len(figures10)} independent panels and 3 layout proofs in {OUT}")


if __name__ == "__main__":
    main()
