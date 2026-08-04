#!/usr/bin/env python3
"""Render current manuscript Figure 3 from distributed source tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import pandas as pd


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
    "B_Plasma": "B/Plasma",
    "Stromal_Endothelial": "Stromal/Endo.",
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


def label(axis, letter: str) -> None:
    axis.text(-0.08, 1.04, letter, transform=axis.transAxes, weight="bold", size=10)


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

    # The distributed CSV basenames retain their frozen historical Figure4
    # labels for checksum provenance; the current manuscript maps these four
    # panels to Figure 3 without changing their contents.
    paths = {
        "a": args.data_dir / "Figure4A_bone_marrow_UMAP_display_source_data_20260716_v2.csv",
        "b": args.data_dir / "Figure4B_subject_timepoint_composition_source_data_20260716_v2.csv",
        "c": args.data_dir / "Figure4C_strict_module_compartment_projection_source_data_20260716_v2.csv",
        "d": args.data_dir / "Figure4D_marrow_support_response_source_data_20260716_v2.csv",
    }
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    a, b, c, d = (pd.read_csv(paths[key]) for key in "abcd")

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 7,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(183 / 25.4, 165 / 25.4))
    ax = axes[0, 0]
    for cell_type in ["Unknown", *CELL_ORDER]:
        part = a.loc[a["broad_cell_type"].eq(cell_type)]
        ax.scatter(
            part["UMAP_1"],
            part["UMAP_2"],
            s=0.25,
            c=CELL_COLORS[cell_type],
            linewidths=0,
            alpha=0.75,
            rasterized=True,
            label=CELL_LABELS[cell_type],
        )
    ax.set_title("Integrated bone-marrow single-cell atlas", loc="left", weight="bold")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(markerscale=5, ncol=3, fontsize=5.2, loc="lower center", bbox_to_anchor=(0.5, -0.33))
    label(ax, "A")

    ax = axes[0, 1]
    med = (
        b.loc[b["broad_cell_type"].isin(CELL_ORDER)]
        .drop_duplicates(["group", "broad_cell_type"])
        .pivot(index="broad_cell_type", columns="group", values="group_median_proportion")
        .reindex(index=CELL_ORDER, columns=GROUP_ORDER)
        * 100
    )
    image = ax.imshow(med, aspect="auto", cmap="Blues", vmin=0)
    for row in range(med.shape[0]):
        for col in range(med.shape[1]):
            ax.text(col, row, f"{med.iloc[row, col]:.1f}", ha="center", va="center", size=5.5)
    ax.set_xticks(range(len(GROUP_ORDER)), [GROUP_LABELS[x] for x in GROUP_ORDER])
    ax.set_yticks(range(len(CELL_ORDER)), [CELL_LABELS[x] for x in CELL_ORDER])
    ax.set_title("Median marrow-cell fractions by clinical group", loc="left", weight="bold")
    fig.colorbar(image, ax=ax, fraction=0.04, pad=0.02, label="%")
    label(ax, "B")

    ax = axes[1, 0]
    base = c.loc[c["variant_id"].eq("hub_top50_mean_all")].copy()
    matrix = (
        base.pivot(
            index="module_color",
            columns="coarse_cell_type",
            values="relative_score_within_module",
        )
        .reindex(columns=CELL_ORDER)
        .sort_index()
    )
    image = ax.imshow(matrix, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(CELL_ORDER)), [CELL_LABELS[x] for x in CELL_ORDER], rotation=45, ha="right")
    ax.set_yticks(range(len(matrix.index)), matrix.index)
    ax.set_title("CD34-derived module projection", loc="left", weight="bold")
    fig.colorbar(image, ax=ax, fraction=0.04, pad=0.02, label="Within-module relative score")
    label(ax, "C")

    ax = axes[1, 1]
    genes = list(dict.fromkeys(d["gene"]))
    matrix = (
        d.pivot(index="gene", columns="coarse_cell_type", values="mean_log1pCP10K")
        .reindex(index=genes, columns=CELL_ORDER)
    )
    cmap = LinearSegmentedColormap.from_list("response", ["#F8F5FA", "#6A1B9A"])
    vmax = float(np.nanmax(matrix.to_numpy()))
    image = ax.imshow(matrix, aspect="auto", cmap=cmap, norm=Normalize(0, vmax))
    detected = (
        d.pivot(index="gene", columns="coarse_cell_type", values="pct_cells_detected")
        .reindex(index=genes, columns=CELL_ORDER)
    )
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            size = 3 + 24 * np.sqrt(max(float(detected.iloc[row, col]), 0) / 100)
            ax.scatter(col, row, s=size, facecolors="none", edgecolors="#333333", linewidths=0.35)
    ax.set_xticks(range(len(CELL_ORDER)), [CELL_LABELS[x] for x in CELL_ORDER], rotation=45, ha="right")
    ax.set_yticks(range(len(genes)), genes)
    ax.set_title("Selected hematopoietic response-gene expression", loc="left", weight="bold")
    fig.colorbar(image, ax=ax, fraction=0.04, pad=0.02, label="Mean log1p(CP10K)")
    label(ax, "D")

    for axis in axes.flat:
        axis.spines[["top", "right"]].set_visible(False)
    fig.subplots_adjust(left=0.10, right=0.95, top=0.95, bottom=0.14, hspace=0.50, wspace=0.38)
    stem = args.output_dir / "Figure_3_bone_marrow_atlas_module_projection"
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)
    print(f"PASS: rendered {stem}")


if __name__ == "__main__":
    main()
