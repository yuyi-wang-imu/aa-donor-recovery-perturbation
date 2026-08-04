#!/usr/bin/env python3
"""Render manuscript Figure 3 with submission-safe gene and merged-category typography.

The distributed source tables retain their frozen historical Figure4 basenames
for provenance. This renderer maps those four source panels to current
manuscript Figure 3 without recalculation or data changes.
"""

from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from PIL import Image


CELL_ORDER = [
    "HSPC",
    "Erythroid",
    "Megakaryocyte",
    "Myeloid",
    "T_NK",
    "B_Plasma",
    "Stromal_Endothelial",
]
CELL_LABELS = {
    "HSPC": "HSPC",
    "Erythroid": "Erythroid",
    "Megakaryocyte": "Megakaryocyte",
    "Myeloid": "Myeloid",
    "T_NK": "T/NK",
    "B_Plasma": "B/plasma",
    "Stromal_Endothelial": "stromal/endothelial",
    "Unknown": "Low-confidence",
}
CELL_COLORS = {
    "HSPC": "#AA3377",
    "Erythroid": "#CC6677",
    "Megakaryocyte": "#EECC66",
    "Myeloid": "#228833",
    "T_NK": "#4477AA",
    "B_Plasma": "#66CCEE",
    "Stromal_Endothelial": "#888888",
    "Unknown": "#D9D9D9",
}
GROUP_ORDER = ["HD", "SAA_baseline", "SAA_3M", "SAA_6M"]
GROUP_LABELS = {"HD": "HD", "SAA_baseline": "Baseline", "SAA_3M": "3M", "SAA_6M": "6M"}
RESPONSE_GENES = ["MPL", "JAK2", "STAT5A", "STAT5B", "PIM1", "BCL2L1"]


def panel_label(axis, letter: str) -> None:
    axis.text(
        -0.12,
        1.04,
        letter,
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        fontweight="bold",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    if not args.data_dir.is_dir():
        raise FileNotFoundError(args.data_dir)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError("Output directory must be new or empty")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    paths = {
        "a": args.data_dir / "Figure4A_bone_marrow_UMAP_display_source_data_20260716_v2.csv",
        "b": args.data_dir / "Figure4B_subject_timepoint_composition_source_data_20260716_v2.csv",
        "c": args.data_dir / "Figure4C_strict_module_compartment_projection_source_data_20260716_v2.csv",
        "d": args.data_dir / "Figure4D_marrow_support_response_source_data_20260716_v2.csv",
        "manifest": args.data_dir / "Figure4_bone_marrow_atlas_projection_source_manifest_20260716_v2.json",
    }
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    a, b, c, d = (pd.read_csv(paths[key]) for key in "abcd")
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
            "font.size": 6.6,
            "axes.linewidth": 0.65,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    fig = plt.figure(figsize=(183 / 25.4, 200 / 25.4), facecolor="white")
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.16, 1.0],
        height_ratios=[1.02, 0.98],
        hspace=0.46,
        wspace=0.38,
    )

    # A: display-only sampling of the frozen RPCA UMAP coordinates.
    ax_a = fig.add_subplot(gs[0, 0])
    draw_order = [
        "Unknown",
        "Stromal_Endothelial",
        "Megakaryocyte",
        "HSPC",
        "B_Plasma",
        "Myeloid",
        "Erythroid",
        "T_NK",
    ]
    for cell_type in draw_order:
        part = a.loc[a["broad_cell_type"].eq(cell_type)]
        ax_a.scatter(
            part["UMAP_1"],
            part["UMAP_2"],
            s=1.15 if cell_type in {"HSPC", "Megakaryocyte", "Stromal_Endothelial"} else 0.65,
            c=CELL_COLORS[cell_type],
            alpha=0.55 if cell_type != "Unknown" else 0.38,
            linewidths=0,
            rasterized=True,
            label=CELL_LABELS[cell_type],
        )
    ax_a.set_xlim(a["UMAP_1"].quantile(0.001), a["UMAP_1"].quantile(0.999))
    ax_a.set_ylim(a["UMAP_2"].quantile(0.001), a["UMAP_2"].quantile(0.999))
    ax_a.set_xticks([])
    ax_a.set_yticks([])
    ax_a.set_xlabel("UMAP 1 (a.u.)", labelpad=2)
    ax_a.set_ylabel("UMAP 2 (a.u.)", labelpad=2)
    ax_a.set_title("Integrated marrow single-cell atlas", loc="left", fontsize=8.2, fontweight="bold", pad=4)
    n_cells = int(manifest.get("frozen_counts", {}).get("cells", 768_617))
    ax_a.text(
        0.01,
        -0.16,
        f"{n_cells:,} cells analyzed; deterministic downsampling for display only",
        transform=ax_a.transAxes,
        fontsize=5.6,
        color="#5F6B76",
        ha="left",
    )
    handles, labels = ax_a.get_legend_handles_labels()
    ax_a.legend(
        handles[::-1],
        labels[::-1],
        ncol=4,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.20),
        fontsize=5.15,
        handletextpad=0.20,
        columnspacing=0.45,
        markerscale=3.8,
    )
    panel_label(ax_a, "A")

    # B: subject-by-timepoint median compartment composition.
    ax_b = fig.add_subplot(gs[0, 1])
    med = (
        b.loc[b["broad_cell_type"].isin(CELL_ORDER)]
        .drop_duplicates(["group", "broad_cell_type"])
        .pivot(index="broad_cell_type", columns="group", values="group_median_proportion")
        .reindex(index=CELL_ORDER, columns=GROUP_ORDER)
        * 100
    )
    comp_cmap = LinearSegmentedColormap.from_list(
        "comp", ["#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"]
    )
    im_b = ax_b.imshow(
        med.to_numpy(),
        aspect="auto",
        cmap=comp_cmap,
        vmin=0,
        vmax=max(60, float(np.nanmax(med.to_numpy()))),
    )
    for row in range(med.shape[0]):
        for col in range(med.shape[1]):
            value = med.iat[row, col]
            ax_b.text(
                col,
                row,
                f"{value:.1f}",
                ha="center",
                va="center",
                fontsize=5.6,
                color="white" if value > 32 else "#17212B",
            )
    ax_b.set_yticks(range(len(CELL_ORDER)))
    ax_b.set_yticklabels([CELL_LABELS[x] for x in CELL_ORDER], fontsize=5.6)
    group_units = (
        b.drop_duplicates(["group", "n_subject_timepoints"])
        .set_index("group")
        .reindex(GROUP_ORDER)["n_subject_timepoints"]
    )
    ax_b.set_xticks(range(len(GROUP_ORDER)))
    ax_b.set_xticklabels(
        [f"{GROUP_LABELS[g]}\nn={int(group_units.loc[g])}" for g in GROUP_ORDER],
        fontsize=5.8,
    )
    ax_b.tick_params(length=0)
    for spine in ax_b.spines.values():
        spine.set_visible(False)
    ax_b.set_title(
        "Median marrow-cell fractions by clinical group",
        loc="left",
        fontsize=8.2,
        fontweight="bold",
        pad=4,
    )
    ax_b.text(
        0.0,
        -0.23,
        "Values are percentages; low-confidence cells excluded from the denominator",
        transform=ax_b.transAxes,
        fontsize=5.45,
        color="#5F6B76",
        ha="left",
    )
    cbar_b = fig.colorbar(im_b, ax=ax_b, fraction=0.045, pad=0.025)
    cbar_b.set_label("Median cells (%)", fontsize=5.8)
    cbar_b.ax.tick_params(labelsize=5.3, length=2)
    panel_label(ax_b, "B")

    # C: strict module projection and release decisions.
    ax_c = fig.add_subplot(gs[1, 0])
    base = c.loc[c["variant_id"].eq("hub_top50_mean_all")].copy()
    modules = sorted(base["module_color"].unique())
    c_pivot = (
        base.pivot(
            index="module_color",
            columns="coarse_cell_type",
            values="relative_score_within_module",
        )
        .reindex(index=modules, columns=CELL_ORDER)
    )
    proj_cmap = LinearSegmentedColormap.from_list(
        "proj", ["#F7FBFF", "#BDD7E7", "#6BAED6", "#2171B5", "#08306B"]
    )
    im_c = ax_c.imshow(c_pivot.to_numpy(), aspect="auto", cmap=proj_cmap, vmin=0, vmax=1)
    ax_c.set_yticks(range(len(modules)))
    ax_c.set_yticklabels(modules, fontsize=5.7)
    ax_c.set_xticks(range(len(CELL_ORDER)))
    ax_c.set_xticklabels(
        [CELL_LABELS[x] for x in CELL_ORDER],
        rotation=45,
        ha="right",
        fontsize=5.2,
    )
    ax_c.tick_params(length=0)
    for spine in ax_c.spines.values():
        spine.set_visible(False)
    decision_index = base.drop_duplicates("module_color").set_index("module_color")
    status_colors = {
        "include_descriptive_context_edge": "#168A5B",
        "downgrade_subject_timepoint_top1_discordance": "#D97706",
        "exclude_instability_below_8_of_9": "#B33A3A",
    }
    status_labels = {
        "include_descriptive_context_edge": "retained",
        "downgrade_subject_timepoint_top1_discordance": "ST discord.",
        "exclude_instability_below_8_of_9": "<8 of 9",
    }
    for row, module in enumerate(modules):
        record = decision_index.loc[module]
        top_cell = record["baseline_top1_cell_type"]
        col = CELL_ORDER.index(top_cell)
        decision = record["edge_decision"]
        ax_c.add_patch(
            Rectangle(
                (col - 0.48, row - 0.48),
                0.96,
                0.96,
                fill=False,
                ec=status_colors[decision],
                lw=1.45,
            )
        )
        text_status = (
            f"{int(record['n_variants_agree_top1'])} of 9  "
            f"{status_labels[decision]}"
        )
        if decision == "downgrade_subject_timepoint_top1_discordance":
            arrow = CELL_LABELS[record["subject_timepoint_top1_cell_type"]]
            text_status += f" \u2192 {arrow}"
        ax_c.text(
            7.18,
            row,
            text_status,
            va="center",
            ha="left",
            fontsize=5.25,
            color=status_colors[decision],
        )
    ax_c.set_xlim(-0.5, 10.8)
    ax_c.set_title("CD34-derived module projection", loc="left", fontsize=8.2, fontweight="bold", pad=4)
    cax_c = ax_c.inset_axes([0.0, -0.31, 0.92, 0.045])
    cbar_c = fig.colorbar(im_c, cax=cax_c, orientation="horizontal")
    cbar_c.set_label("Within-module relative mean score", fontsize=5.8)
    cbar_c.ax.tick_params(labelsize=5.2, length=2)
    panel_label(ax_c, "C")

    # D: marrow-side gene-expression context. Human gene symbols are italic.
    ax_d = fig.add_subplot(gs[1, 1])
    d = d.copy()
    d["gene"] = pd.Categorical(d["gene"], categories=RESPONSE_GENES, ordered=True)
    d["coarse_cell_type"] = pd.Categorical(
        d["coarse_cell_type"], categories=CELL_ORDER, ordered=True
    )
    d = d.sort_values(["gene", "coarse_cell_type"])
    x_map = {cell_type: i for i, cell_type in enumerate(CELL_ORDER)}
    y_map = {gene: i for i, gene in enumerate(RESPONSE_GENES)}
    x = d["coarse_cell_type"].map(x_map).astype(float).to_numpy()
    y = d["gene"].map(y_map).astype(float).to_numpy()
    mean_expr = d["mean_log1pCP10K"].to_numpy()
    detection = d["pct_cells_detected"].to_numpy()
    norm_d = Normalize(vmin=0, vmax=float(np.nanmax(mean_expr)))
    dot_cmap = LinearSegmentedColormap.from_list(
        "dot", ["#F4F1F7", "#C2A5CF", "#7B3294", "#3F007D"]
    )
    sizes = 8 + detection * 3.8
    sc = ax_d.scatter(
        x,
        y,
        s=sizes,
        c=mean_expr,
        cmap=dot_cmap,
        norm=norm_d,
        ec="#5B4A66",
        lw=0.35,
    )
    ax_d.set_xlim(-0.6, len(CELL_ORDER) - 0.4)
    ax_d.set_ylim(len(RESPONSE_GENES) - 0.45, -0.55)
    ax_d.set_xticks(range(len(CELL_ORDER)))
    ax_d.set_xticklabels(
        [CELL_LABELS[x] for x in CELL_ORDER],
        rotation=45,
        ha="right",
        fontsize=5.2,
    )
    ax_d.set_yticks(range(len(RESPONSE_GENES)))
    ax_d.set_yticklabels(RESPONSE_GENES, fontsize=6.0, fontstyle="italic")
    ax_d.grid(color="#E7E7E7", lw=0.45)
    ax_d.set_axisbelow(True)
    ax_d.tick_params(length=0)
    for spine in ax_d.spines.values():
        spine.set_visible(False)
    ax_d.set_title("Marrow-side response-gene context", loc="left", fontsize=8.2, fontweight="bold", pad=4)
    cax_d = ax_d.inset_axes([0.0, -0.31, 0.61, 0.045])
    cbar_d = fig.colorbar(sc, cax=cax_d, orientation="horizontal")
    cbar_d.set_label("Mean log1p(CP10K)", fontsize=5.8)
    cbar_d.ax.tick_params(labelsize=5.2, length=2)
    size_handles = [
        ax_d.scatter([], [], s=8 + pct * 3.8, c="white", ec="#5B4A66", lw=0.5)
        for pct in (5, 20, 40)
    ]
    ax_d.legend(
        size_handles,
        ["5%", "20%", "40%"],
        title="Cells detected",
        ncol=3,
        loc="upper left",
        bbox_to_anchor=(0.65, -0.265),
        fontsize=5.0,
        title_fontsize=5.2,
        handletextpad=0.15,
        columnspacing=0.35,
        borderaxespad=0,
    )
    panel_label(ax_d, "D")

    fig.subplots_adjust(left=0.09, right=0.97, top=0.965, bottom=0.15)
    stem = args.output_dir / "Figure_3_bone_marrow_atlas_gene_format_fixed_v2_20260726"
    fig.savefig(stem.with_suffix(".svg"), format="svg", dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), format="pdf", dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".png"), format="png", dpi=600, facecolor="white")
    plt.close(fig)

    # Write through in-memory buffers so non-ASCII Windows paths remain safe.
    with Image.open(stem.with_suffix(".png")) as image:
        buffer = io.BytesIO()
        image.convert("RGB").save(
            buffer, format="PNG", dpi=(600, 600), optimize=True
        )
        png_payload = buffer.getvalue()
    stem.with_suffix(".png").write_bytes(png_payload)
    with Image.open(stem.with_suffix(".png")) as image:
        buffer = io.BytesIO()
        image.convert("RGB").save(
            buffer,
            format="TIFF",
            dpi=(600, 600),
            compression="tiff_lzw",
        )
        tiff_payload = buffer.getvalue()
    stem.with_suffix(".tif").write_bytes(tiff_payload)
    print(f"PASS: rendered {stem}")


if __name__ == "__main__":
    main()
