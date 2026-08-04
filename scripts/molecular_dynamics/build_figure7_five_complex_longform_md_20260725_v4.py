#!/usr/bin/env python3
"""Layout-corrected v4 of the five-complex 100-ns long-form MD figures.

Relative to v3, this script changes only typography and spacing:
- the legend occupies a dedicated band above panels A-B;
- long column titles are shortened without changing metric definitions;
- compact complex labels prevent clipping in the final-20-ns table.
No data, axes, smoothing, time range, or summary calculation is changed.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np


VERSION = "20260725_v4"
BASE_SCRIPT = pathlib.Path(__file__).with_name(
    "build_figure7_five_complex_longform_md_20260725_v3.py"
)


def load_base():
    spec = importlib.util.spec_from_file_location("figure7_v3_base", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import plotting base: {BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BASE = load_base()
CASES = BASE.CASES
DISPLAY = BASE.DISPLAY
SHORT_DISPLAY = BASE.SHORT_DISPLAY
COLORS = BASE.COLORS
TIME_SERIES = BASE.TIME_SERIES
RMSF_SERIES = BASE.RMSF_SERIES

SHORT_TITLE = {
    "protein_backbone_rmsd": "Protein-backbone RMSD",
    "ligand_heavy_atom_rmsd_after_protein_fit": "Ligand RMSD after protein fit",
    "global_protein_ligand_minimum_distance": "Minimum protein-ligand distance",
    "ligand_heavy_atom_proximity_fraction_lt_0p45nm": (
        "Ligand proximity fraction (<0.45 nm)"
    ),
}


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def make_overview(frame, limits, output_dir):
    figure = plt.figure(
        figsize=(183 / 25.4, 238 / 25.4),
        facecolor="white",
    )
    outer = figure.add_gridspec(
        4,
        1,
        height_ratios=[1.0, 1.0, 0.78, 0.72],
        left=0.105,
        right=0.985,
        bottom=0.055,
        top=0.865,
        hspace=0.72,
    )
    top = outer[:2].subgridspec(2, 2, wspace=0.28, hspace=0.48)
    overlay_axes = [
        figure.add_subplot(top[row, column])
        for row in range(2)
        for column in range(2)
    ]
    for letter, axis, metric in zip("ABCD", overlay_axes, TIME_SERIES):
        BASE.plot_overlay(axis, frame, metric, limits)
        axis.set_title(
            SHORT_TITLE[metric],
            loc="left",
            fontweight="bold",
            pad=4,
        )
        BASE.add_panel_letter(axis, letter)

    handles, labels = overlay_axes[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.55, 0.91),
        ncol=3,
        columnspacing=1.35,
        handlelength=2.8,
    )

    rmsf_grid = outer[2].subgridspec(1, 5, wspace=0.30)
    rmsf_axes = [figure.add_subplot(rmsf_grid[0, index]) for index in range(5)]
    for index, (axis, case_id) in enumerate(zip(rmsf_axes, CASES)):
        BASE.style_axis(axis)
        x, y = BASE.series(frame, case_id, RMSF_SERIES)
        axis.plot(x, y, color=COLORS[case_id], linewidth=0.35, alpha=0.20)
        axis.plot(
            x,
            BASE.BASE.rolling_median(y, 7),
            color=COLORS[case_id],
            linewidth=1.05,
        )
        axis.set_ylim(*limits[RMSF_SERIES])
        axis.set_title(
            SHORT_DISPLAY[case_id],
            fontsize=6.2,
            color=COLORS[case_id],
            fontweight="bold",
            pad=3,
        )
        axis.set_xlabel("Residue index")
        if index == 0:
            axis.set_ylabel("C-alpha RMSF (nm)")
            BASE.add_panel_letter(axis, "E")
        else:
            axis.set_ylabel("")

    summary_axis = figure.add_subplot(outer[3])
    normalized, labels_text = BASE.final20_table(frame)
    summary_axis.imshow(
        normalized,
        aspect="auto",
        cmap=plt.matplotlib.colors.LinearSegmentedColormap.from_list(
            "neutral_summary", ["#F6F7F7", "#D4E3EB", "#739CB5"]
        ),
        vmin=0,
        vmax=1,
    )
    summary_axis.set_xticks(range(4))
    summary_axis.set_xticklabels(
        [
            "Protein RMSD\nmean (nm)",
            "Ligand RMSD\nmean (nm)",
            "Minimum distance\nmedian (nm)",
            "Proximity fraction\nmean",
        ]
    )
    summary_axis.set_yticks(range(5))
    summary_axis.set_yticklabels([SHORT_DISPLAY[case_id] for case_id in CASES])
    for row in range(5):
        for column in range(4):
            summary_axis.text(
                column,
                row,
                labels_text[row][column],
                ha="center",
                va="center",
                fontsize=6.7,
                fontweight="bold",
                color=("#FFFFFF" if normalized[row, column] > 0.62 else "#26343C"),
            )
    summary_axis.set_title(
        "Final 20 ns quantitative summary",
        loc="left",
        fontweight="bold",
        pad=5,
    )
    summary_axis.tick_params(length=0, pad=3)
    for spine in summary_axis.spines.values():
        spine.set_visible(False)
    BASE.add_panel_letter(summary_axis, "F")
    summary_axis.text(
        1.0,
        -0.25,
        "Cell shading is scaled within each column; printed values are the reported quantities.",
        transform=summary_axis.transAxes,
        ha="right",
        va="top",
        fontsize=5.3,
        color="#65717A",
    )

    figure.text(
        0.105,
        0.976,
        "Standardized 100-ns molecular-dynamics assessment of five representative complexes",
        fontsize=10.5,
        fontweight="bold",
        va="center",
    )
    figure.text(
        0.105,
        0.952,
        (
            "All time-series panels span the complete 0-100 ns trajectories; "
            "grey shading denotes the final 20 ns."
        ),
        fontsize=6.2,
        color="#56616A",
    )
    stem = f"Figure_7_five_complexes_100ns_MD_overview_{VERSION}"
    files = BASE.save_all(figure, output_dir, stem)
    plt.close(figure)
    return files


def make_full_matrix(frame, limits, output_dir):
    figure = plt.figure(
        figsize=(183 / 25.4, 230 / 25.4),
        facecolor="white",
    )
    grid = figure.add_gridspec(
        5,
        4,
        left=0.14,
        right=0.985,
        bottom=0.055,
        top=0.90,
        wspace=0.38,
        hspace=0.58,
    )
    for row, case_id in enumerate(CASES):
        for column, metric in enumerate(TIME_SERIES):
            axis = figure.add_subplot(grid[row, column])
            BASE.style_axis(axis)
            x, y = BASE.series(frame, case_id, metric)
            axis.plot(x, y, color=COLORS[case_id], linewidth=0.25, alpha=0.12)
            axis.plot(
                x,
                BASE.BASE.rolling_median(y, 201),
                color=COLORS[case_id],
                linewidth=1.15,
            )
            axis.axvspan(80, 100, color="#AEB9C3", alpha=0.14, linewidth=0)
            axis.set_xlim(0, 100)
            axis.set_xticks([0, 25, 50, 75, 100])
            axis.set_ylim(*limits[metric])
            axis.set_xlabel("Time (ns)")
            if row == 0:
                axis.set_title(
                    SHORT_TITLE[metric],
                    fontsize=6.7,
                    fontweight="bold",
                    pad=4,
                )
            if column == 0:
                axis.set_ylabel(BASE.PANEL_YLABEL[metric])
                axis.text(
                    -0.52,
                    0.5,
                    DISPLAY[case_id],
                    transform=axis.transAxes,
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=7.1,
                    fontweight="bold",
                    color=COLORS[case_id],
                )
                axis.text(
                    -0.52,
                    1.10,
                    chr(ord("A") + row),
                    transform=axis.transAxes,
                    fontsize=10,
                    fontweight="bold",
                    color="#111111",
                )
            else:
                axis.set_ylabel("")
    figure.text(
        0.14,
        0.972,
        "Complete 0-100 ns trajectories displayed separately for all five complexes",
        fontsize=10.5,
        fontweight="bold",
        va="center",
    )
    figure.text(
        0.14,
        0.944,
        (
            "Thin traces show frame-level values; thick traces show centered rolling medians; "
            "grey shading denotes 80-100 ns."
        ),
        fontsize=6.1,
        color="#56616A",
    )
    stem = f"Figure_7_five_complexes_100ns_MD_full_matrix_{VERSION}"
    files = BASE.save_all(figure, output_dir, stem)
    plt.close(figure)
    return files


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-csv", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError(
            f"Refusing to write into existing directory: {args.output_dir}"
        )
    if not args.source_csv.is_file():
        raise FileNotFoundError(args.source_csv)

    BASE.configure_style()
    frame = BASE.BASE.load_and_validate(args.source_csv)
    limits = BASE.BASE.global_limits(frame)
    args.output_dir.mkdir(parents=True, exist_ok=False)

    overview = make_overview(frame, limits, args.output_dir)
    full_matrix = make_full_matrix(frame, limits, args.output_dir)
    outputs = [*overview.values(), *full_matrix.values()]
    manifest = {
        "version": VERSION,
        "change_from_v3": (
            "Typography and spacing only; no change to source data, axes, "
            "time range, smoothing, or final-20-ns calculations."
        ),
        "source_csv": str(args.source_csv),
        "source_csv_sha256": sha256(args.source_csv),
        "cases": [
            {"case_id": case_id, "display": DISPLAY[case_id], "duration_ns": 100}
            for case_id in CASES
        ],
        "outputs": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in outputs
        },
    }
    manifest_path = args.output_dir / (
        f"Figure_7_five_complex_longform_manifest_{VERSION}.json"
    )
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
