#!/usr/bin/env python3
"""Static release checks that do not require licensed or large external data."""

from __future__ import annotations

import argparse
import csv
import io
import re
from pathlib import Path, PurePosixPath
from zipfile import ZipFile

from build_manifest import (
    OUTPUT as MANIFEST_PATH,
    ROOT,
    canonical_bytes,
    digest_bytes,
    manifest_rows,
    render_manifest,
)


MANIFEST_FIELDS = ["relative_path", "bytes", "sha256"]


def fail(message: str) -> None:
    raise SystemExit(message)


def parse_manifest(data: bytes) -> list[dict[str, str]]:
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        fail(f"Manifest is not UTF-8: {exc}")
    reader = csv.DictReader(io.StringIO(text, newline=""), delimiter="\t")
    if reader.fieldnames != MANIFEST_FIELDS:
        fail(f"Invalid manifest header: {reader.fieldnames}")
    rows = list(reader)
    paths = [row["relative_path"] for row in rows]
    if len(paths) != len(set(paths)):
        fail("Manifest contains duplicate relative paths")
    return rows


def verify_archive(archive_path: Path) -> None:
    with ZipFile(archive_path, "r") as archive:
        files = [info for info in archive.infolist() if not info.is_dir()]
        manifest_entries = [
            info
            for info in files
            if PurePosixPath(info.filename).name == "MANIFEST.tsv"
        ]
        if len(manifest_entries) != 1:
            fail(
                "Release archive must contain exactly one MANIFEST.tsv; "
                f"found {len(manifest_entries)}"
            )
        manifest_entry = manifest_entries[0]
        prefix = manifest_entry.filename[: -len("MANIFEST.tsv")]
        rows = parse_manifest(archive.read(manifest_entry))
        expected_names = {prefix + row["relative_path"] for row in rows}
        expected_names.add(manifest_entry.filename)
        actual_names = {info.filename for info in files}
        if actual_names != expected_names:
            missing = sorted(expected_names - actual_names)
            extra = sorted(actual_names - expected_names)
            fail(f"Archive coverage mismatch: missing={missing}, extra={extra}")
        for row in rows:
            name = prefix + row["relative_path"]
            data = archive.read(name)
            if len(data) != int(row["bytes"]) or digest_bytes(data) != row["sha256"]:
                fail(f"Archive manifest mismatch: {row['relative_path']}")
    print(f"PASS: archive manifest verified for {len(rows)} files")


def validate_repository() -> int:
    required = [
        ".gitattributes",
        ".gitignore",
        "LICENSE",
        "README.md",
        "CITATION.cff",
        "DATA_AND_LICENSES.md",
        "FIGURE_SOURCE_MAP.tsv",
        "GITHUB_ZENODO_RELEASE_CHECKLIST.md",
        "WORKFLOW_ORDER.tsv",
        "config/analysis_parameters.tsv",
        "config/random_seeds.tsv",
        "environment/requirements.txt",
        "environment/geneformer_gpu_environment_20260804_v1.yml",
        "environment/install_r_packages.R",
        "environment/r_packages.tsv",
        "scripts/run_workflow.py",
        "scripts/build_manifest.py",
        "scripts/build_release_zip.py",
        "scripts/validate_r_parse.R",
        "scripts/geneformer/70_bmc_geneformer_prepare_balanced_cd34_public_20260804_v1.R",
        "scripts/geneformer/71_bmc_geneformer_donor_mvp_20260802_v1.py",
        "scripts/geneformer/85_bmc_geneformer_overexpression_mvp_20260802_v3.py",
        "scripts/geneformer/93_bmc_geneformer_bidirectional_audit_20260802_v2.py",
        "scripts/geneformer/109_bmc_geneformer_state_specific_perturbation_20260802_v1.py",
        "scripts/geneformer_figures/build_bmc_geneformer_panels_20260803_v1.py",
        "scripts/molecular_dynamics/build_figure7_five_complex_standard_md_20260725_v2.py",
        "scripts/molecular_dynamics/build_figure7_five_complex_longform_md_20260725_v3.py",
        "scripts/molecular_dynamics/build_figure7_five_complex_longform_md_20260725_v4.py",
        "scripts/transcriptomics_wgcna/build_figure2_panel_label_pvalue_geneitalic_exactcontent_20260726_v7.R",
        "scripts/figure_packaging/build_figure4_5_layout_fixes_20260726_v1.py",
        "scripts/molecular_dynamics/build_figure7_five_complex_longform_md_20260726_v6.py",
        "scripts/transcriptomics_wgcna/build_figure2_candidate_pool_CD34_context_typography_20260726_v7.py",
        "scripts/figure_packaging/assemble_supplementary_figures_s1_s8_20260727_v2.py",
        "derived_data/molecular_dynamics/Figure_7_five_system_md_source_data_20260725_v5_pbc_audited.csv",
        "derived_data/bone_marrow/Figure4A_bone_marrow_UMAP_display_source_data_20260716_v2.csv",
        "derived_data/bone_marrow/Figure4B_subject_timepoint_composition_source_data_20260716_v2.csv",
        "derived_data/bone_marrow/Figure4C_strict_module_compartment_projection_source_data_20260716_v2.csv",
        "derived_data/bone_marrow/Figure4D_marrow_support_response_source_data_20260716_v2.csv",
        "derived_data/computational_perturbation/Figure8_standard_panelB_response_landscape_source_20260725_v4.csv",
        "derived_data/computational_perturbation/Figure8_standard_panelC_ranked_response_source_20260725_v4.csv",
        "derived_data/computational_perturbation/Figure8_standard_panelD_response_breadth_source_20260725_v4.csv",
        "derived_data/computational_perturbation/FigureS11_matched_control_calibration_source_20260725_v4.csv",
        "derived_data/computational_perturbation/FigureS11_matched_control_covariate_balance_source_20260725_v4.csv",
        "derived_data/computational_perturbation/FigureS11_matched_control_pooled_null_source_20260725_v4.csv",
    ]
    missing = [name for name in required if not (ROOT / name).is_file()]
    if missing:
        fail(f"Missing release files: {', '.join(missing)}")

    with (ROOT / "WORKFLOW_ORDER.tsv").open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        steps = list(csv.DictReader(handle, delimiter="\t"))
    missing_entries = [
        row["canonical_entry"]
        for row in steps
        if not (ROOT / row["canonical_entry"]).is_file()
    ]
    if missing_entries:
        fail(f"Missing canonical entries: {', '.join(missing_entries)}")

    seed_text = (ROOT / "config/random_seeds.tsv").read_text(encoding="utf-8")
    if "sctenifold_formal\t20260724" not in seed_text:
        fail("Formal scTenifoldKnk common seed is not registered")
    if "geneformer_perturbation\t20260802" not in seed_text:
        fail("Geneformer common seed is not registered")

    forbidden_suffixes = {
        ".xtc", ".trr", ".tpr", ".edr", ".cpt", ".dcd", ".nc",
        ".safetensors", ".pt", ".pth", ".ckpt", ".pem", ".key", ".p12", ".pfx",
    }
    forbidden_files = [
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file() and path.suffix.lower() in forbidden_suffixes
    ]
    if forbidden_files:
        fail(f"Forbidden large/sensitive release files: {forbidden_files}")

    high_confidence_secret_patterns = [
        re.compile(r"AKIA[0-9A-Z]{16}"),
        re.compile(r"gh[pousr]_[A-Za-z0-9]{30,}"),
        re.compile(r"-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----"),
    ]
    secret_hits = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in {".py", ".r", ".md", ".tsv", ".csv", ".json", ".yml", ".yaml", ".cff"}:
            continue
        text = path.read_text(encoding="utf-8-sig", errors="ignore")
        if any(pattern.search(text) for pattern in high_confidence_secret_patterns):
            secret_hits.append(path.relative_to(ROOT).as_posix())
    if secret_hits:
        fail(f"Possible secrets/private keys found: {secret_hits}")

    for metadata_name in ("README.md", "CITATION.cff"):
        metadata = (ROOT / metadata_name).read_text(encoding="utf-8-sig")
        historical_doi = "10.5281/" + "zenodo.21644540"
        historical_repo = "aa-prescription-" + "hematopoietic-workflow"
        if historical_doi in metadata or historical_repo in metadata:
            fail(f"Historical repository/DOI leaked into new release metadata: {metadata_name}")

    if not MANIFEST_PATH.is_file():
        fail("Missing release file: MANIFEST.tsv")
    rows = parse_manifest(canonical_bytes(MANIFEST_PATH))
    expected_rows = manifest_rows()
    expected_paths = {row[0] for row in expected_rows}
    actual_paths = {row["relative_path"] for row in rows}
    if actual_paths != expected_paths:
        missing = sorted(expected_paths - actual_paths)
        extra = sorted(actual_paths - expected_paths)
        fail(f"Manifest coverage mismatch: missing={missing}, extra={extra}")
    for row in rows:
        path = ROOT / row["relative_path"]
        data = canonical_bytes(path)
        if len(data) != int(row["bytes"]) or digest_bytes(data) != row["sha256"]:
            fail(f"Manifest mismatch: {row['relative_path']}")
    if canonical_bytes(MANIFEST_PATH) != render_manifest(expected_rows):
        fail("MANIFEST.tsv is not in canonical path order or is stale")

    print(
        f"PASS: {len(steps)} canonical entries, {len(rows)} manifest files, "
        "and release metadata verified"
    )
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--archive",
        type=Path,
        help="also verify an existing release ZIP against its embedded manifest",
    )
    args = parser.parse_args()
    validate_repository()
    if args.archive is not None:
        verify_archive(args.archive.resolve())


if __name__ == "__main__":
    main()
