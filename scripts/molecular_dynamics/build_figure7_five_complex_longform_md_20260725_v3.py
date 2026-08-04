#!/usr/bin/env python3
"""Build long-axis, standardized 100-ns MD figures for five complexes.

This plotting-only script reads the frozen PBC-audited source table used by the
20260725 v2 review bundle.  It does not change trajectories, measurements,
manuscript text, or any previous figure.  Two complementary review figures are
written:

1. a submission-oriented overview with four long time-series panels, five
   complex-specific C-alpha RMSF panels, and a final-20-ns numeric summary;
2. a complete five-row time-series matrix in which every complex is displayed
   separately against the full 0-100 ns interval.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import pathlib

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


VERSION = "20260725_v3"
TAIL_START_NS = 80.0
BASE_SCRIPT = pathlib.Path(__file__).with_name(
    "build_figure7_five_complex_standard_md_20260725_v2.py"
)


def load_base():
    spec = importlib.util.spec_from_file_location("figure7_v2_base", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import plotting base: {BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BASE = load_base()
CASES = BASE.CASES
DISPLAY = BASE.DISPLAY
COLORS = BASE.COLORS
TIME_SERIES = BASE.TIME_SERIES
RMSF_SERIES = BASE.RMSF_SERIES

SHORT_DISPLAY = {
    "PARP1_sesamin": "PARP1-sesamin",
    "KIT_3O_Methylorobol": "KIT-3'-O-Me",
    "CDK6_3O_Methylorobol": "CDK6-3'-O-Me",
    "SYK_isofucosterol": "SYK-isofucosterol",
    "GSK3B_linarin": "GSK3B-linarin",
}

PANEL_TITLE = {
    "protein_backbone_rmsd": "Protein-backbone RMSD",
    "ligand_heavy_atom_rmsd_after_protein_fit": (
        "Ligand RMSD after protein-backbone fit"
    ),
    "global_protein_ligand_minimum_distance": (
        "Global protein-ligand minimum distance"
    ),
    "ligand_heavy_atom_proximity_fraction_lt_0p45nm": (
        "Ligand-heavy-atom proximity fraction (<0.45 nm)"
    ),
}

PANEL_YLABEL = {
    "protein_backbone_rmsd": "RMSD (nm)",
    "ligand_heavy_atom_rmsd_after_protein_fit": "RMSD (nm)",
    "global_protein_ligand_minimum_distance": "Distance (nm)",
    "ligand_heavy_atom_proximity_fraction_lt_0p45nm": "Fraction",
}


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def configure_style() -> None:
    BASE.configure_style()
    mpl.rcParams.update(
        {
            "font.size": 7.0,
            "axes.titlesize": 8.0,
            "axes.labelsize": 7.2,
            "xtick.labelsize": 6.4,
            "ytick.labelsize": 6.4,
            "legend.fontsize": 6.5,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
        }
    )


def series(
    frame: pd.DataFrame, case_id: str, metric: str
) -> tuple[np.ndarray, np.ndarray]:
    data = BASE.values(frame, case_id, metric)
    return (
        data["x"].to_numpy(dtype=float),
        data["value"].to_numpy(dtype=float),
    )


def style_axis(axis: plt.Axes) -> None:
    axis.grid(color="#DDE5EA", linewidth=0.42, alpha=0.82)
    axis.set_axisbelow(True)
    axis.tick_params(pad=2.0)


def add_panel_letter(axis: plt.Axes, letter: str) -> None:
    axis.text(
        -0.12,
        1.08,
        letter,
        transform=axis.transAxes,
        fontsize=10.5,
        fontweight="bold",
        va="bottom",
        ha="left",
        color="#111111",
    )


def save_all(
    figure: plt.Figure, output_dir: pathlib.Path, stem: str
) -> dict[str, pathlib.Path]:
    paths = {
        "svg": output_dir / f"{stem}.svg",
        "pdf": output_dir / f"{stem}.pdf",
        "png": output_dir / f"{stem}.png",
        "tiff": output_dir / f"{stem}.tiff",
    }
    for path in paths.values():
        if path.exists():
            raise FileExistsError(path)
    figure.savefig(paths["svg"], facecolor="white")
    figure.savefig(paths["pdf"], facecolor="white")
    figure.savefig(paths["png"], dpi=600, facecolor="white")
    figure.savefig(
        paths["tiff"],
        dpi=600,
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    return paths


def plot_overlay(
    axis: plt.Axes,
    frame: pd.DataFrame,
    metric: str,
    limits: dict[str, tuple[float, float]],
) -> None:
    style_axis(axis)
    for case_id in CASES:
        x, y = series(frame, case_id, metric)
        color = COLORS[case_id]
        axis.plot(x, y, color=color, linewidth=0.24, alpha=0.08)
        axis.plot(
            x,
            BASE.rolling_median(y, 201),
            color=color,
            linewidth=1.25,
            label=DISPLAY[case_id],
        )
    axis.axvspan(80, 100, color="#AEB9C3", alpha=0.14, linewidth=0)
    axis.set_xlim(0, 100)
    axis.set_xticks([0, 20, 40, 60, 80, 100])
    axis.set_ylim(*limits[metric])
    axis.set_xlabel("Time (ns)")
    axis.set_ylabel(PANEL_YLABEL[metric])
    axis.set_title(PANEL_TITLE[metric], loc="left", fontweight="bold", pad=4)


def final20_table(frame: pd.DataFrame) -> tuple[np.ndarray, list[list[str]]]:
    metric_order = (
        "protein_backbone_rmsd",
        "ligand_heavy_atom_rmsd_after_protein_fit",
        "global_protein_ligand_minimum_distance",
        "ligand_heavy_atom_proximity_fraction_lt_0p45nm",
    )
    raw = []
    text = []
    for case_id in CASES:
        stats = BASE.final20_stats(frame, case_id)
        row = [stats[metric] for metric in metric_order]
        raw.append(row)
        text.append(
            [
                f"{row[0]:.3f}",
                f"{row[1]:.3f}",
                f"{row[2]:.3f}",
                f"{row[3]:.3f}",
            ]
        )
    array = np.asarray(raw, dtype=float)
    normalized = np.zeros_like(array)
    for column in range(array.shape[1]):
        values = array[:, column]
        span = float(values.max() - values.min())
        normalized[:, column] = (
            0.5 if span == 0 else (values - values.min()) / span
        )
    return normalized, text


def make_overview(
    frame: pd.DataFrame,
    limits: dict[str, tuple[float, float]],
    output_dir: pathlib.Path,
) -> dict[str, pathlib.Path]:
    figure = plt.figure(
        figsize=(183 / 25.4, 235 / 25.4),
        facecolor="white",
    )
    outer = figure.add_gridspec(
        4,
        1,
        height_ratios=[1.0, 1.0, 0.78, 0.72],
        left=0.09,
        right=0.985,
        bottom=0.055,
        top=0.91,
        hspace=0.72,
    )
    top = outer[:2].subgridspec(2, 2, wspace=0.28, hspace=0.48)
    overlay_axes = [
        figure.add_subplot(top[row, column])
        for row in range(2)
        for column in range(2)
    ]
    for letter, axis, metric in zip(
        "ABCD", overlay_axes, TIME_SERIES
    ):
        plot_overlay(axis, frame, metric, limits)
        add_panel_letter(axis, letter)

    handles, labels = overlay_axes[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.55, 0.932),
        ncol=3,
        columnspacing=1.25,
        handlelength=2.5,
    )

    rmsf_grid = outer[2].subgridspec(1, 5, wspace=0.30)
    rmsf_axes = [figure.add_subplot(rmsf_grid[0, index]) for index in range(5)]
    for index, (axis, case_id) in enumerate(zip(rmsf_axes, CASES)):
        style_axis(axis)
        x, y = series(frame, case_id, RMSF_SERIES)
        axis.plot(x, y, color=COLORS[case_id], linewidth=0.35, alpha=0.20)
        axis.plot(
            x,
            BASE.rolling_median(y, 7),
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
            add_panel_letter(axis, "E")
        else:
            axis.set_ylabel("")

    summary_axis = figure.add_subplot(outer[3])
    normalized, labels_text = final20_table(frame)
    summary_axis.imshow(
        normalized,
        aspect="auto",
        cmap=mpl.colors.LinearSegmentedColormap.from_list(
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
    summary_axis.set_yticklabels([DISPLAY[case_id] for case_id in CASES])
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
    add_panel_letter(summary_axis, "F")
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
        0.09,
        0.975,
        "Standardized 100-ns molecular-dynamics assessment of five representative complexes",
        fontsize=10.5,
        fontweight="bold",
        va="center",
    )
    figure.text(
        0.09,
        0.948,
        (
            "All time-series panels span the complete 0-100 ns trajectories; "
            "grey shading denotes the final 20 ns."
        ),
        fontsize=6.2,
        color="#56616A",
    )
    stem = f"Figure_7_five_complexes_100ns_MD_overview_{VERSION}"
    files = save_all(figure, output_dir, stem)
    plt.close(figure)
    return files


def make_full_matrix(
    frame: pd.DataFrame,
    limits: dict[str, tuple[float, float]],
    output_dir: pathlib.Path,
) -> dict[str, pathlib.Path]:
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
        top=0.91,
        wspace=0.38,
        hspace=0.58,
    )
    for row, case_id in enumerate(CASES):
        for column, metric in enumerate(TIME_SERIES):
            axis = figure.add_subplot(grid[row, column])
            style_axis(axis)
            x, y = series(frame, case_id, metric)
            axis.plot(x, y, color=COLORS[case_id], linewidth=0.25, alpha=0.12)
            axis.plot(
                x,
                BASE.rolling_median(y, 201),
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
                    PANEL_TITLE[metric],
                    fontsize=6.7,
                    fontweight="bold",
                    pad=4,
                )
            if column == 0:
                axis.set_ylabel(PANEL_YLABEL[metric])
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
    files = save_all(figure, output_dir, stem)
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

    configure_style()
    frame = BASE.load_and_validate(args.source_csv)
    limits = BASE.global_limits(frame)
    args.output_dir.mkdir(parents=True, exist_ok=False)

    overview = make_overview(frame, limits, args.output_dir)
    full_matrix = make_full_matrix(frame, limits, args.output_dir)

    outputs = [*overview.values(), *full_matrix.values()]
    manifest = {
        "version": VERSION,
        "purpose": (
            "Long-axis review redesign using only the frozen PBC-audited "
            "five-system 100-ns source table."
        ),
        "source_csv": str(args.source_csv),
        "source_csv_sha256": sha256(args.source_csv),
        "cases": [
            {"case_id": case_id, "display": DISPLAY[case_id], "duration_ns": 100}
            for case_id in CASES
        ],
        "metrics": list(TIME_SERIES) + [RMSF_SERIES],
        "display": {
            "overview": (
                "four complete 0-100 ns overlay panels, five separate C-alpha "
                "RMSF panels, and final-20-ns numeric summary"
            ),
            "full_matrix": (
                "five separate complex rows by four complete 0-100 ns time-series columns"
            ),
            "final_window_ns": [80, 100],
        },
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
