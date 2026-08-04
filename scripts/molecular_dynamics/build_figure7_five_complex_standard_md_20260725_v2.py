#!/usr/bin/env python3
"""Build standardized 100-ns MD plates for all five representative complexes.

The script is plotting-only.  It reads the frozen PBC-audited five-system
source table, validates the case/metric inventory and duration, and writes a
new review bundle without changing the manuscript or any existing figure.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import pathlib

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


VERSION = "20260725_v2"
TAIL_START_NS = 80.0
EXPECTED_END_NS = 100.0
TIME_SERIES = (
    "protein_backbone_rmsd",
    "ligand_heavy_atom_rmsd_after_protein_fit",
    "global_protein_ligand_minimum_distance",
    "ligand_heavy_atom_proximity_fraction_lt_0p45nm",
)
RMSF_SERIES = "ca_rmsf_by_residue"
CASES = (
    "PARP1_sesamin",
    "KIT_3O_Methylorobol",
    "CDK6_3O_Methylorobol",
    "SYK_isofucosterol",
    "GSK3B_linarin",
)
DISPLAY = {
    "PARP1_sesamin": "PARP1-sesamin",
    "KIT_3O_Methylorobol": "KIT-3'-O-methylorobol",
    "CDK6_3O_Methylorobol": "CDK6-3'-O-methylorobol",
    "SYK_isofucosterol": "SYK-isofucosterol",
    "GSK3B_linarin": "GSK3B-linarin",
}
COLORS = {
    "PARP1_sesamin": "#3977B7",
    "KIT_3O_Methylorobol": "#2F9D83",
    "CDK6_3O_Methylorobol": "#D9792B",
    "SYK_isofucosterol": "#7B63B5",
    "GSK3B_linarin": "#7C9A38",
}
METRIC_TITLE = {
    "protein_backbone_rmsd": "Protein-backbone RMSD",
    "ligand_heavy_atom_rmsd_after_protein_fit": (
        "Ligand RMSD after protein fit"
    ),
    "global_protein_ligand_minimum_distance": (
        "Protein-ligand minimum distance"
    ),
    "ligand_heavy_atom_proximity_fraction_lt_0p45nm": (
        "Ligand-heavy-atom proximity (<0.45 nm)"
    ),
    RMSF_SERIES: "C-alpha RMSF",
}
MASTER_METRIC_TITLE = {
    "protein_backbone_rmsd": "Protein RMSD",
    "ligand_heavy_atom_rmsd_after_protein_fit": (
        "Ligand RMSD\n(after protein fit)"
    ),
    "global_protein_ligand_minimum_distance": "Minimum distance",
    "ligand_heavy_atom_proximity_fraction_lt_0p45nm": (
        "Proximity fraction\n(<0.45 nm)"
    ),
    RMSF_SERIES: "C-alpha RMSF",
}
METRIC_YLABEL = {
    "protein_backbone_rmsd": "RMSD (nm)",
    "ligand_heavy_atom_rmsd_after_protein_fit": "RMSD (nm)",
    "global_protein_ligand_minimum_distance": "Distance (nm)",
    "ligand_heavy_atom_proximity_fraction_lt_0p45nm": "Fraction",
    RMSF_SERIES: "RMSF (nm)",
}


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 6.5,
            "axes.titlesize": 7.0,
            "axes.labelsize": 6.2,
            "xtick.labelsize": 5.5,
            "ytick.labelsize": 5.5,
            "axes.linewidth": 0.55,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.major.size": 2.4,
            "ytick.major.size": 2.4,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def load_and_validate(source_csv: pathlib.Path) -> pd.DataFrame:
    frame = pd.read_csv(source_csv)
    required_columns = {"case_id", "series", "x", "value", "x_unit", "value_unit"}
    missing = required_columns.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    if set(frame["case_id"].unique()) != set(CASES):
        raise ValueError("The source table does not contain exactly the five cases.")
    required_series = set(TIME_SERIES) | {RMSF_SERIES}
    for case_id in CASES:
        case = frame.loc[frame["case_id"] == case_id]
        if set(case["series"].unique()) != required_series:
            raise ValueError(f"Metric inventory mismatch: {case_id}")
        for metric in TIME_SERIES:
            values = case.loc[case["series"] == metric].sort_values("x")
            if len(values) != 10001:
                raise ValueError(f"{case_id}/{metric}: expected 10001 frames")
            x = values["x"].to_numpy(dtype=float)
            y = values["value"].to_numpy(dtype=float)
            if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
                raise ValueError(f"{case_id}/{metric}: non-finite values")
            if not np.all(np.diff(x) > 0):
                raise ValueError(f"{case_id}/{metric}: non-increasing time")
            if not np.isclose(x[0], 0.0, atol=0.02):
                raise ValueError(f"{case_id}/{metric}: does not start at 0 ns")
            if not np.isclose(x[-1], EXPECTED_END_NS, atol=0.02):
                raise ValueError(f"{case_id}/{metric}: does not end at 100 ns")
        rmsf = case.loc[case["series"] == RMSF_SERIES, "value"].to_numpy(float)
        if rmsf.size == 0 or not np.all(np.isfinite(rmsf)) or np.any(rmsf < 0):
            raise ValueError(f"{case_id}: invalid RMSF")
    if not np.all(np.isfinite(frame[["x", "value"]].to_numpy(float))):
        raise ValueError("Non-finite values in source table.")
    return frame


def values(frame: pd.DataFrame, case_id: str, metric: str) -> pd.DataFrame:
    result = frame.loc[
        (frame["case_id"] == case_id) & (frame["series"] == metric),
        ["x", "value"],
    ].sort_values("x")
    if result.empty:
        raise ValueError(f"Empty series: {case_id}/{metric}")
    return result


def rolling_median(y: np.ndarray, window: int) -> np.ndarray:
    return (
        pd.Series(y)
        .rolling(window=window, center=True, min_periods=max(3, window // 5))
        .median()
        .to_numpy()
    )


def final20_stats(frame: pd.DataFrame, case_id: str) -> dict[str, float]:
    stats: dict[str, float] = {}
    for metric in TIME_SERIES:
        data = values(frame, case_id, metric)
        tail = data.loc[data["x"] >= TAIL_START_NS, "value"].to_numpy(float)
        stats[metric] = (
            float(np.median(tail))
            if metric == "global_protein_ligand_minimum_distance"
            else float(np.mean(tail))
        )
    rmsf = values(frame, case_id, RMSF_SERIES)["value"].to_numpy(float)
    stats["median_ca_rmsf"] = float(np.median(rmsf))
    return stats


def global_limits(frame: pd.DataFrame) -> dict[str, tuple[float, float]]:
    limits: dict[str, tuple[float, float]] = {}
    for metric in (
        "protein_backbone_rmsd",
        "ligand_heavy_atom_rmsd_after_protein_fit",
        RMSF_SERIES,
    ):
        array = frame.loc[frame["series"] == metric, "value"].to_numpy(float)
        upper = float(np.nanpercentile(array, 99.7)) * 1.08
        limits[metric] = (0.0, upper)
    distance = frame.loc[
        frame["series"] == "global_protein_ligand_minimum_distance", "value"
    ].to_numpy(float)
    padding = max(0.005, (np.nanpercentile(distance, 99.7) -
                          np.nanpercentile(distance, 0.3)) * 0.08)
    limits["global_protein_ligand_minimum_distance"] = (
        float(np.nanpercentile(distance, 0.3) - padding),
        float(np.nanpercentile(distance, 99.7) + padding),
    )
    limits["ligand_heavy_atom_proximity_fraction_lt_0p45nm"] = (0.0, 1.02)
    return limits


def style_axis(axis: plt.Axes) -> None:
    axis.grid(color="#DEE5EA", linewidth=0.38, alpha=0.82)
    axis.set_axisbelow(True)


def plot_metric(
    axis: plt.Axes,
    frame: pd.DataFrame,
    case_id: str,
    metric: str,
    limits: dict[str, tuple[float, float]],
    show_xlabel: bool = True,
) -> None:
    data = values(frame, case_id, metric)
    x = data["x"].to_numpy(float)
    y = data["value"].to_numpy(float)
    color = COLORS[case_id]
    style_axis(axis)
    if metric == RMSF_SERIES:
        axis.plot(x, y, color=color, linewidth=0.45, alpha=0.28)
        axis.plot(x, rolling_median(y, 7), color=color, linewidth=1.0)
        if show_xlabel:
            axis.set_xlabel("Residue index")
    else:
        axis.plot(x, y, color=color, linewidth=0.35, alpha=0.18)
        axis.plot(x, rolling_median(y, 201), color=color, linewidth=1.0)
        axis.axvspan(80, 100, color="#AEB9C3", alpha=0.16, linewidth=0)
        axis.set_xlim(0, 100)
        if show_xlabel:
            axis.set_xlabel("Time (ns)")
    axis.set_ylim(*limits[metric])
    axis.set_ylabel(METRIC_YLABEL[metric])


def save_figure(
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
    figure.savefig(paths["png"], dpi=300, facecolor="white")
    figure.savefig(
        paths["tiff"],
        dpi=600,
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    return paths


def make_case_plate(
    frame: pd.DataFrame,
    case_id: str,
    limits: dict[str, tuple[float, float]],
    output_dir: pathlib.Path,
) -> dict[str, pathlib.Path]:
    figure = plt.figure(
        figsize=(183 / 25.4, 118 / 25.4),
        facecolor="white",
    )
    grid = figure.add_gridspec(
        2,
        3,
        left=0.085,
        right=0.985,
        bottom=0.105,
        top=0.86,
        wspace=0.42,
        hspace=0.58,
    )
    axes = [figure.add_subplot(grid[row, col]) for row in range(2) for col in range(3)]
    metrics = (
        "protein_backbone_rmsd",
        "ligand_heavy_atom_rmsd_after_protein_fit",
        "global_protein_ligand_minimum_distance",
        "ligand_heavy_atom_proximity_fraction_lt_0p45nm",
        RMSF_SERIES,
    )
    for index, (axis, metric) in enumerate(zip(axes, metrics)):
        plot_metric(axis, frame, case_id, metric, limits)
        axis.set_title(METRIC_TITLE[metric], loc="left", fontweight="bold", pad=3)
        axis.text(
            -0.16,
            1.08,
            chr(ord("A") + index),
            transform=axis.transAxes,
            fontsize=9,
            fontweight="bold",
            va="bottom",
        )

    summary = axes[5]
    summary.axis("off")
    summary.text(
        0.0,
        1.02,
        "F",
        transform=summary.transAxes,
        fontsize=9,
        fontweight="bold",
        va="bottom",
    )
    summary.text(
        0.08,
        0.95,
        "Final 20 ns summary",
        fontsize=7.0,
        fontweight="bold",
        va="top",
    )
    stats = final20_stats(frame, case_id)
    lines = (
        ("Protein RMSD, mean", stats["protein_backbone_rmsd"], "nm"),
        (
            "Ligand RMSD, mean",
            stats["ligand_heavy_atom_rmsd_after_protein_fit"],
            "nm",
        ),
        (
            "Minimum distance, median",
            stats["global_protein_ligand_minimum_distance"],
            "nm",
        ),
        (
            "Proximity fraction, mean",
            stats["ligand_heavy_atom_proximity_fraction_lt_0p45nm"],
            "",
        ),
        ("C-alpha RMSF, median", stats["median_ca_rmsf"], "nm"),
    )
    y = 0.78
    for label, value, unit in lines:
        summary.text(0.08, y, label, fontsize=5.9, color="#39434A")
        summary.text(
            0.95,
            y,
            f"{value:.3f}{(' ' + unit) if unit else ''}",
            fontsize=6.2,
            fontweight="bold",
            ha="right",
            color=COLORS[case_id],
        )
        y -= 0.135
    summary.text(
        0.08,
        0.05,
        "Thin lines: frame-level values; thick lines: centered rolling medians.",
        fontsize=5.0,
        color="#66727A",
        wrap=True,
    )

    figure.text(
        0.085,
        0.945,
        f"{DISPLAY[case_id]}: standardized 100-ns molecular-dynamics descriptors",
        fontsize=9.5,
        fontweight="bold",
        va="center",
    )
    stem = f"{case_id}_standardized_100ns_MD_{VERSION}"
    files = save_figure(figure, output_dir, stem)
    plt.close(figure)
    return files


def make_master_plate(
    frame: pd.DataFrame,
    limits: dict[str, tuple[float, float]],
    output_dir: pathlib.Path,
) -> dict[str, pathlib.Path]:
    figure = plt.figure(
        figsize=(183 / 25.4, 225 / 25.4),
        facecolor="white",
    )
    grid = figure.add_gridspec(
        5,
        5,
        left=0.13,
        right=0.99,
        bottom=0.055,
        top=0.91,
        wspace=0.48,
        hspace=0.55,
    )
    metrics = (
        "protein_backbone_rmsd",
        "ligand_heavy_atom_rmsd_after_protein_fit",
        "global_protein_ligand_minimum_distance",
        "ligand_heavy_atom_proximity_fraction_lt_0p45nm",
        RMSF_SERIES,
    )
    for row, case_id in enumerate(CASES):
        for col, metric in enumerate(metrics):
            axis = figure.add_subplot(grid[row, col])
            plot_metric(
                axis,
                frame,
                case_id,
                metric,
                limits,
                show_xlabel=(row == len(CASES) - 1),
            )
            if row == 0:
                axis.set_title(
                    MASTER_METRIC_TITLE[metric],
                    fontsize=5.9,
                    fontweight="bold",
                    pad=4,
                )
            if col != 0:
                axis.set_ylabel("")
            if row != len(CASES) - 1:
                axis.set_xticklabels([])
            if col == 0:
                axis.text(
                    -0.5,
                    0.5,
                    DISPLAY[case_id],
                    transform=axis.transAxes,
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=6.4,
                    fontweight="bold",
                    color=COLORS[case_id],
                )
                axis.text(
                    -0.5,
                    1.17,
                    chr(ord("A") + row),
                    transform=axis.transAxes,
                    fontsize=9,
                    fontweight="bold",
                    va="top",
                    ha="center",
                    color="#111111",
                )
    figure.text(
        0.13,
        0.965,
        "Standardized 100-ns molecular-dynamics descriptors for five representative complexes",
        fontsize=9.5,
        fontweight="bold",
        va="center",
    )
    figure.text(
        0.13,
        0.935,
        "Each row is one complex; all rows use the same metric definitions and column-specific scales.",
        fontsize=5.8,
        color="#56616A",
    )
    stem = f"Figure_7_five_complexes_standardized_100ns_MD_{VERSION}"
    files = save_figure(figure, output_dir, stem)
    plt.close(figure)
    return files


def write_summary_csv(
    frame: pd.DataFrame, output_dir: pathlib.Path
) -> pathlib.Path:
    path = output_dir / f"Figure_7_five_complex_final20_summary_{VERSION}.csv"
    if path.exists():
        raise FileExistsError(path)
    fields = (
        "case_id",
        "complex",
        "protein_backbone_rmsd_mean_nm",
        "ligand_rmsd_after_protein_fit_mean_nm",
        "global_minimum_distance_median_nm",
        "proximity_fraction_mean",
        "ca_rmsf_median_nm",
    )
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for case_id in CASES:
            stats = final20_stats(frame, case_id)
            writer.writerow(
                {
                    "case_id": case_id,
                    "complex": DISPLAY[case_id],
                    "protein_backbone_rmsd_mean_nm": (
                        f"{stats['protein_backbone_rmsd']:.9g}"
                    ),
                    "ligand_rmsd_after_protein_fit_mean_nm": (
                        f"{stats['ligand_heavy_atom_rmsd_after_protein_fit']:.9g}"
                    ),
                    "global_minimum_distance_median_nm": (
                        f"{stats['global_protein_ligand_minimum_distance']:.9g}"
                    ),
                    "proximity_fraction_mean": (
                        f"{stats['ligand_heavy_atom_proximity_fraction_lt_0p45nm']:.9g}"
                    ),
                    "ca_rmsf_median_nm": f"{stats['median_ca_rmsf']:.9g}",
                }
            )
    return path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-csv", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"Refusing existing output directory: {args.output_dir}")
    if not args.source_csv.is_file():
        raise FileNotFoundError(args.source_csv)

    configure_style()
    frame = load_and_validate(args.source_csv)
    limits = global_limits(frame)
    args.output_dir.mkdir(parents=True, exist_ok=False)

    master_files = make_master_plate(frame, limits, args.output_dir)
    case_files = {
        case_id: make_case_plate(frame, case_id, limits, args.output_dir)
        for case_id in CASES
    }
    summary_csv = write_summary_csv(frame, args.output_dir)

    outputs = [*master_files.values(), summary_csv]
    for files in case_files.values():
        outputs.extend(files.values())
    manifest = {
        "version": VERSION,
        "source_csv": str(args.source_csv),
        "source_csv_sha256": sha256(args.source_csv),
        "cases": [
            {"case_id": case_id, "display": DISPLAY[case_id], "duration_ns": 100}
            for case_id in CASES
        ],
        "metric_definitions": {
            "protein_backbone_rmsd": (
                "Protein-backbone RMSD from the PBC-corrected trajectory."
            ),
            "ligand_heavy_atom_rmsd_after_protein_fit": (
                "Ligand-heavy-atom RMSD after least-squares protein-backbone fit."
            ),
            "global_protein_ligand_minimum_distance": (
                "Global minimum heavy-atom distance between protein and ligand."
            ),
            "ligand_heavy_atom_proximity_fraction_lt_0p45nm": (
                "Fraction of ligand heavy atoms within 0.45 nm of any protein heavy atom."
            ),
            "ca_rmsf_by_residue": "C-alpha RMSF by protein-specific residue index.",
        },
        "display": {
            "master_plate": "five rows by five standardized metric columns",
            "per_complex_plate": "six panels per complex, including a final-20-ns summary",
            "final_window_ns": [80, 100],
        },
        "outputs": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in outputs
        },
    }
    manifest_path = args.output_dir / f"Figure_7_standardized_MD_manifest_{VERSION}.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
