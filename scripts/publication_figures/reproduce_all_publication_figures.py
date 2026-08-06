#!/usr/bin/env python3
"""Rebuild all manuscript figures and verify them against frozen references.

This is the publication-replay layer. It uses only repository-distributed
derived tables and approved author-generated publication intermediates. It
does not claim to rerun licensed database exports, molecular docking, raw MD
trajectories, or Geneformer model inference.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from PIL import Image

from assemble_static_publication_figures import align_panel_letters, default_font


def run(command: list[str], cwd: Path, env: dict[str, str] | None = None) -> None:
    print("RUN:", subprocess.list2cmdline(command), flush=True)
    subprocess.run(command, cwd=cwd, env=env, check=True)


def copy(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    if destination.exists():
        raise FileExistsError(destination)
    shutil.copy2(source, destination)


def erase_top_title_band(path: Path, rows: int) -> None:
    """Erase only a frozen top title band while preserving all lower pixels."""

    with Image.open(path) as image:
        dpi = image.info.get("dpi", (600, 600))
        original = np.asarray(image.convert("RGB"), dtype=np.uint8).copy()
    modified = original.copy()
    modified[:rows, :, :] = 255
    temporary = path.with_name(path.stem + ".title-band.tmp.png")
    Image.fromarray(modified, mode="RGB").save(temporary, dpi=dpi)
    with Image.open(temporary) as image:
        verified = np.asarray(image.convert("RGB"), dtype=np.uint8)
    if verified.shape != original.shape:
        raise RuntimeError(f"Title-band postprocess changed canvas dimensions: {path}")
    if not np.array_equal(verified[rows:, :, :], original[rows:, :, :]):
        raise RuntimeError(f"Title-band postprocess changed panel pixels: {path}")
    if not np.all(verified[:rows, :, :] == 255):
        raise RuntimeError(f"Title-band postprocess did not produce a pure-white band: {path}")
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument(
        "--reference-root",
        type=Path,
        help="Defaults to repository/reference_outputs.",
    )
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[2]
    output_root = args.output_root.resolve()
    if output_root.exists():
        raise FileExistsError(f"Refusing to overwrite existing output root: {output_root}")
    main_dir = output_root / "main_figures"
    supplementary_dir = output_root / "supplementary_figures"
    work_dir = output_root / "work"
    main_dir.mkdir(parents=True)
    supplementary_dir.mkdir()
    work_dir.mkdir()

    static_dir = work_dir / "static"
    run(
        [
            sys.executable,
            "-B",
            str(repo / "scripts/publication_figures/assemble_static_publication_figures.py"),
            "--asset-dir",
            str(repo / "derived_data/publication_intermediates"),
            "--output-dir",
            str(static_dir),
        ],
        repo,
    )
    for number in (1, 2, 4, 5, 6):
        copy(static_dir / f"Figure_{number}.png", main_dir / f"Figure_{number}.png")
    erase_top_title_band(main_dir / "Figure_5.png", 300)

    figure3_raw = work_dir / "figure3_raw"
    run(
        [
            sys.executable,
            "-B",
            str(repo / "scripts/bone_marrow_single_cell/render_bone_marrow_figure3_gene_format_fixed_v5.py"),
            str(repo / "derived_data/bone_marrow"),
            str(figure3_raw),
        ],
        repo,
    )
    align_panel_letters(
        figure3_raw / "Figure_3_bone_marrow_atlas_gene_format_fixed_v2_20260726.png",
        main_dir / "Figure_3.png",
        3,
        default_font(),
    )

    figure7_dir = work_dir / "figure7"
    run(
        [
            sys.executable,
            "-B",
            str(repo / "scripts/molecular_dynamics/build_figure7_five_complex_longform_md_20260726_v6.py"),
            "--source-csv",
            str(
                repo
                / "derived_data/molecular_dynamics/Figure_7_five_system_md_source_data_20260725_v5_pbc_audited.csv"
            ),
            "--output-dir",
            str(figure7_dir),
            "--publication-final",
        ],
        repo,
    )
    copy(
        figure7_dir
        / "Figure_7_five_complexes_100ns_MD_overview_20260806_no_top_title_final.png",
        main_dir / "Figure_7.png",
    )

    figure8_dir = work_dir / "figure8"
    figure8_dir.mkdir()
    run(
        [
            args.rscript,
            str(
                repo
                / "scripts/computational_perturbation/build_figure8_standard_and_s11_sensitivity_candidate_v5.R"
            ),
            str(repo / "derived_data/computational_perturbation"),
            str(figure8_dir),
        ],
        repo,
    )
    copy(
        figure8_dir
        / "Figure_8_standard_scTenifoldKnk_network_responses_pvalue_format_20260726_v12.png",
        main_dir / "Figure_8.png",
    )
    copy(
        figure8_dir
        / "Supplementary_Figure_S11_matched_control_sensitivity_pvalue_format_20260726_v6.png",
        supplementary_dir / "Figure_S8.png",
    )

    geneformer_dir = work_dir / "geneformer"
    geneformer_env = os.environ.copy()
    geneformer_env["AA_GENEFORMER_DATA_DIR"] = str(repo / "derived_data/geneformer")
    geneformer_env["AA_GENEFORMER_FIGURE_OUT"] = str(geneformer_dir)
    run(
        [
            sys.executable,
            "-B",
            str(repo / "scripts/geneformer_figures/build_bmc_geneformer_panels_20260803_v6.py"),
        ],
        repo,
        geneformer_env,
    )
    run(
        [
            sys.executable,
            "-B",
            str(repo / "scripts/geneformer_figures/convert_bmc_geneformer_panels_rgb_20260803_v7.py"),
        ],
        repo,
        geneformer_env,
    )
    copy(
        geneformer_dir / "Figure_9_Layout_Proof_20260803_v7.png",
        main_dir / "Figure_9.png",
    )
    copy(
        geneformer_dir / "Figure_S9_Layout_Proof_20260803_v7.png",
        supplementary_dir / "Figure_S9.png",
    )
    copy(
        geneformer_dir / "Figure_S10_Layout_Proof_20260803_v7.png",
        supplementary_dir / "Figure_S10.png",
    )

    supplementary_package_dir = work_dir / "supplementary_package"
    run(
        [
            sys.executable,
            "-B",
            str(
                repo
                / "scripts/figure_packaging/assemble_supplementary_figures_s1_s8_20260727_v2.py"
            ),
            str(
                repo
                / "derived_data/publication_intermediates/supplementary_pages"
            ),
            str(supplementary_package_dir),
        ],
        repo,
    )
    copy(
        supplementary_package_dir / "Supplementary_Figures_S1-S8.pdf",
        supplementary_dir / "Supplementary_Figures_S1-S8.pdf",
    )

    reference_root = (args.reference_root or repo / "reference_outputs").resolve()
    report = output_root / "publication_figure_verification.json"
    run(
        [
            sys.executable,
            "-B",
            str(repo / "scripts/publication_figures/verify_publication_figures.py"),
            "--actual-root",
            str(output_root),
            "--reference-root",
            str(reference_root),
            "--report",
            str(report),
        ],
        repo,
    )
    print(f"PASS: publication replay verified in {output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
