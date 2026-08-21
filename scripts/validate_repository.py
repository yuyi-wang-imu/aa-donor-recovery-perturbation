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
        "REPRODUCIBILITY_MATRIX.tsv",
        "PUBLICATION_ASSET_CHECKSUMS.tsv",
        "CURRENT_MANUSCRIPT_FIGURE_SOURCE_MAP.tsv",
        "CURRENT_MANUSCRIPT_REPRODUCIBILITY_MATRIX.tsv",
        "CURRENT_MANUSCRIPT_ASSET_CHECKSUMS.tsv",
        "CURRENT_MANUSCRIPT_SUPPLEMENTARY_FIGURE_MAP.tsv",
        "GITHUB_ZENODO_RELEASE_CHECKLIST.md",
        "WORKFLOW_ORDER.tsv",
        "config/analysis_parameters.tsv",
        "config/input_paths.example.tsv",
        "config/random_seeds.tsv",
        "environment/README.md",
        "environment/requirements.txt",
        "environment/publication_replay_python_20260806.txt",
        "environment/figure7_exact_renderer_python_20260806.txt",
        "environment/publication_replay_r_packages_20260806.tsv",
        "environment/geneformer_gpu_environment_20260804_v1.yml",
        "environment/install_r_packages.R",
        "environment/r_packages.tsv",
        "scripts/run_workflow.py",
        "scripts/build_manifest.py",
        "scripts/build_release_zip.py",
        "scripts/validate_r_parse.R",
        "scripts/geneformer/70_aa_geneformer_prepare_balanced_cd34_public_20260804_v1.R",
        "scripts/geneformer/71_aa_geneformer_donor_mvp_20260802_v1.py",
        "scripts/geneformer/85_aa_geneformer_overexpression_mvp_20260802_v3.py",
        "scripts/geneformer/93_aa_geneformer_bidirectional_audit_20260802_v2.py",
        "scripts/geneformer/109_aa_geneformer_state_specific_perturbation_20260802_v1.py",
        "scripts/geneformer_figures/build_publication_geneformer_panels_20260803_v6.py",
        "scripts/geneformer_figures/convert_publication_geneformer_panels_rgb_20260803_v7.py",
        "scripts/molecular_dynamics/build_figure7_five_complex_standard_md_20260725_v2.py",
        "scripts/molecular_dynamics/build_figure7_five_complex_longform_md_20260725_v3.py",
        "scripts/molecular_dynamics/build_figure7_five_complex_longform_md_20260725_v4.py",
        "scripts/transcriptomics_wgcna/build_figure2_panel_label_pvalue_geneitalic_exactcontent_20260726_v7.R",
        "scripts/figure_packaging/build_figure4_5_layout_fixes_20260726_v1.py",
        "scripts/molecular_dynamics/build_figure7_five_complex_longform_md_20260726_v6.py",
        "scripts/molecular_dynamics/render_current_figure8_md.R",
        "scripts/transcriptomics_wgcna/build_figure2_candidate_pool_CD34_context_typography_20260726_v7.py",
        "scripts/figure_packaging/assemble_supplementary_figures_s1_s8_20260727_v2.py",
        "scripts/figure_packaging/build_current_supplementary_visual_corrections.R",
        "scripts/publication_figures/assemble_static_publication_figures.py",
        "scripts/publication_figures/reproduce_all_publication_figures.py",
        "scripts/publication_figures/verify_publication_figures.py",
        "scripts/publication_tables/verify_submission_assets.py",
        "scripts/publication_tables/verify_current_submission_assets.py",
        "derived_data/molecular_dynamics/Figure_7_five_system_md_source_data_20260725_v5_pbc_audited.csv",
        "derived_data/molecular_dynamics/current_figure8/Figure8_five_candidates_time_series_source.tsv.gz",
        "derived_data/molecular_dynamics/current_figure8/Figure8_five_candidates_ca_rmsf_source.tsv",
        "derived_data/molecular_dynamics/current_figure8/Figure8_final20ns_quantitative_summary.tsv",
        "derived_data/molecular_dynamics/current_figure8/README.md",
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
        "derived_data/publication_intermediates/README.md",
        "reference_outputs/README.md",
        "reference_outputs/current_manuscript/README.md",
    ]
    required.extend(
        f"reference_outputs/main_figures/Figure_{number}.png"
        for number in range(1, 10)
    )
    required.extend(
        f"reference_outputs/current_manuscript/Figure_{number}.png"
        for number in range(1, 9)
    )
    required.append("reference_outputs/current_manuscript/Graphical_Abstract.png")
    required.extend(
        f"reference_outputs/supplementary_figures/Figure_{suffix}.png"
        for suffix in ("S8", "S9", "S10")
    )
    required.extend(
        [
            "derived_data/publication_intermediates/supplementary_pages/"
            "Supplementary_Figure_S1_20260726.pdf",
            "derived_data/publication_intermediates/supplementary_pages/"
            "Supplementary_Figure_S2_20260726.pdf",
            "derived_data/publication_intermediates/supplementary_pages/"
            "Supplementary_Figure_S3_DetectionCoverageOnly_20260726.pdf",
            "derived_data/publication_intermediates/supplementary_pages/"
            "Supplementary_Figure_S4_20260726.pdf",
            "derived_data/publication_intermediates/supplementary_pages/"
            "Supplementary_Figure_S5_20260726.pdf",
            "derived_data/publication_intermediates/supplementary_pages/"
            "Supplementary_Figure_S6_ModuleLocalizationRobustness_20260726.pdf",
            "derived_data/publication_intermediates/supplementary_pages/"
            "Supplementary_Figure_S7_DockingMatrixOnly_20260726.pdf",
            "derived_data/publication_intermediates/supplementary_pages/"
            "Supplementary_Figure_S8_MatchedControlSensitivity_20260726.pdf",
        ]
    )
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

    oversized = [
        (path.relative_to(ROOT).as_posix(), path.stat().st_size)
        for path in ROOT.rglob("*")
        if path.is_file()
        and ".git" not in path.relative_to(ROOT).parts
        and path.stat().st_size > 50 * 1024 * 1024
    ]
    if oversized:
        fail(f"Files exceed the 50 MiB release ceiling: {oversized}")

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

    readme = (ROOT / "README.md").read_text(encoding="utf-8-sig")
    required_readme_terms = [
        "https://github.com/yuyi-wang-imu/aa-donor-recovery-perturbation",
        "python3 scripts/validate_repository.py",
        "panel rendering",
        "final RGB conversion",
        "without recomputing analytical results",
        "reproduce_all_publication_figures.py",
        "REPRODUCIBILITY_MATRIX.tsv",
        "verify_submission_assets.py",
        "CURRENT_MANUSCRIPT_FIGURE_SOURCE_MAP.tsv",
        "CURRENT_MANUSCRIPT_REPRODUCIBILITY_MATRIX.tsv",
        "CURRENT_MANUSCRIPT_ASSET_CHECKSUMS.tsv",
        "verify_current_submission_assets.py",
        "render_current_figure8_md.R",
        "build_current_supplementary_visual_corrections.R",
    ]
    missing_readme_terms = [term for term in required_readme_terms if term not in readme]
    if missing_readme_terms:
        fail(f"README is missing publication guidance: {missing_readme_terms}")
    forbidden_readme_terms = [
        "UPLOAD_READINESS_AUDIT_20260804_v1.md",
        "PROPOSED_UPLOAD_MANIFEST_20260804_v1.tsv",
        "remains an author decision",
    ]
    leaked_readme_terms = [term for term in forbidden_readme_terms if term in readme]
    if leaked_readme_terms:
        fail(f"README contains internal or stale guidance: {leaked_readme_terms}")

    golden_hashes = {
        "reference_outputs/main_figures/Figure_1.png": "8FC1DB2D457AE795688AE5449B97985B0CE905049B19D9FA0C845A6C8FDD3B98",
        "reference_outputs/main_figures/Figure_2.png": "7F803D3150BA552A7BE4BF6D966C3D96CCE066ED2BF1BF04D3700B2CA47CFBC1",
        "reference_outputs/main_figures/Figure_3.png": "DB2F65B4A27E914ADF46CA2C8615D327BE59F58B10249AE0FD3B6696E1FEA6D1",
        "reference_outputs/main_figures/Figure_4.png": "117ED322C86BFD888902FF200EA9B5B07534DAAA870CEC3A3FEAC35BB4C8FCF3",
        "reference_outputs/main_figures/Figure_5.png": "922BD60B7C8931ECD23FAB7528A7E423BFA83DE32C69F37840222222E72F9CBE",
        "reference_outputs/main_figures/Figure_6.png": "93E4DBAB25FAA507D32FA9DC377536560BE92D7682723B3365B616952E5DBA1B",
        "reference_outputs/main_figures/Figure_7.png": "BD4215EEB8B14474EF0ED1625D6AD85028776253EF421C828D1319B7893F2A6B",
        "reference_outputs/main_figures/Figure_8.png": "EF41370F994774E0C46C1447E48EDF77D57E0EA4C19C757D96CAED7102A51731",
        "reference_outputs/main_figures/Figure_9.png": "5CB383358B113F04FBE7F4A827817548A3E98534DF27F0EB50081DA0BACCC9CC",
        "reference_outputs/supplementary_figures/Figure_S8.png": "A18A1691FAFADD746C55AEE08316B87126A0F61830A2C445363E91EF4E4DDBD5",
        "reference_outputs/supplementary_figures/Figure_S9.png": "16CC237FDEC974676C3FD3D9ADF1FB309DBB4DDB9E49BB3FC431A57361758E74",
        "reference_outputs/supplementary_figures/Figure_S10.png": "5924F2DDA18953F1647B6BE292D8181347409967291EA6EFF595D9039F57B5DA",
        "derived_data/publication_intermediates/supplementary_pages/Supplementary_Figure_S1_20260726.pdf": "CB1A9FF441576C6F3B1F30E355362518C49FA52E1BC103199116629CE2D32C8A",
        "derived_data/publication_intermediates/supplementary_pages/Supplementary_Figure_S2_20260726.pdf": "4716A73B273823EA39A23B79CB9638D2083293ED1FE12AD79432E6A1F1E7ED3B",
        "derived_data/publication_intermediates/supplementary_pages/Supplementary_Figure_S3_DetectionCoverageOnly_20260726.pdf": "15DD182C7AC7C1E272AC9AA6C426B33E6129AC53706688FA3B2B2BADD42F7EF0",
        "derived_data/publication_intermediates/supplementary_pages/Supplementary_Figure_S4_20260726.pdf": "42ADF2988485037E8926EA0147D5D5C0C2FABE5C814019CECEF106D1B4A0E419",
        "derived_data/publication_intermediates/supplementary_pages/Supplementary_Figure_S5_20260726.pdf": "E06F2D9FC3FD51065EAB194380676DC71C16F087640B4C6BC690DA7A718E5E53",
        "derived_data/publication_intermediates/supplementary_pages/Supplementary_Figure_S6_ModuleLocalizationRobustness_20260726.pdf": "B6500C86F945908F45F8C9687D9DFB2DDA2E194148D20B1ECD4C19794048B7D1",
        "derived_data/publication_intermediates/supplementary_pages/Supplementary_Figure_S7_DockingMatrixOnly_20260726.pdf": "B1DCDDBEEB16E6702447D31C26A56D7266DD1F0F4BA7D03247C785E7956661AE",
        "derived_data/publication_intermediates/supplementary_pages/Supplementary_Figure_S8_MatchedControlSensitivity_20260726.pdf": "737AE4F7C595A85FF03523889A6589C18F68B186A26F8FA8BC006860C13CD3C5",
        "reference_outputs/current_manuscript/Figure_1.png": "7F803D3150BA552A7BE4BF6D966C3D96CCE066ED2BF1BF04D3700B2CA47CFBC1",
        "reference_outputs/current_manuscript/Figure_2.png": "B8166744F62B11B0071A2EC01D6E790C0CB393BA5054AC9945F954BFD83ACC22",
        "reference_outputs/current_manuscript/Figure_3.png": "117ED322C86BFD888902FF200EA9B5B07534DAAA870CEC3A3FEAC35BB4C8FCF3",
        "reference_outputs/current_manuscript/Figure_4.png": "922BD60B7C8931ECD23FAB7528A7E423BFA83DE32C69F37840222222E72F9CBE",
        "reference_outputs/current_manuscript/Figure_5.png": "5CB383358B113F04FBE7F4A827817548A3E98534DF27F0EB50081DA0BACCC9CC",
        "reference_outputs/current_manuscript/Figure_6.png": "EF41370F994774E0C46C1447E48EDF77D57E0EA4C19C757D96CAED7102A51731",
        "reference_outputs/current_manuscript/Figure_7.png": "F620E61BC83AF0C8FA292C21F185D6A23B4BAE02C41F3FAD99331E474043FB39",
        "reference_outputs/current_manuscript/Figure_8.png": "046BAD345A05FFC90D043B742E24A2AD82537CF1C61AD77BCB3A1BC607B6367B",
        "reference_outputs/current_manuscript/Graphical_Abstract.png": "B792B411849F8199031DE6170FCE0CDDBED76717FBD38F6DBF309C9821FB9E07",
        "derived_data/molecular_dynamics/current_figure8/Figure8_five_candidates_time_series_source.tsv.gz": "EFFAF8AA397AF7DDF1365FE6A8246B2B46280E8E532FB175587A5131DF5D5710",
        "derived_data/molecular_dynamics/current_figure8/Figure8_five_candidates_ca_rmsf_source.tsv": "DA6344A272C1C688641BC8F7B4F7BC568A658103BFB80F2F57BB7764FE7E6D13",
        "derived_data/molecular_dynamics/current_figure8/Figure8_final20ns_quantitative_summary.tsv": "63AFC6229EC931481AD209387107475FA06ADC4D4AFE15E09F9D87E23CF0B764",
    }
    golden_mismatches = [
        name
        for name, expected in golden_hashes.items()
        if digest_bytes((ROOT / name).read_bytes()) != expected
    ]
    if golden_mismatches:
        fail(f"Frozen publication references changed: {golden_mismatches}")

    with (ROOT / "PUBLICATION_ASSET_CHECKSUMS.tsv").open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        publication_assets = list(csv.DictReader(handle, delimiter="\t"))
    expected_submission_paths = {
        *(f"02_Main_Figures/Figure_{index}.png" for index in range(1, 10)),
        "04_Additional_Files/Additional_file_1_Supplementary_Tables_S1-S6_PublicationReady.xlsx",
        "04_Additional_Files/Additional_file_2_Supplementary_Figures_S1-S10.pdf",
        "04_Additional_Files/Additional_file_3_Supplementary_Table_S7.xlsx",
        "04_Additional_Files/Additional_file_4_Supplementary_Table_S8.xlsx",
        "04_Additional_Files/Additional_file_5_Supplementary_Table_S10.xlsx",
        "04_Additional_Files/Additional_file_6_Supplementary_Table_S9.xlsx",
    }
    actual_submission_paths = {
        row.get("relative_submission_path", "") for row in publication_assets
    }
    if len(publication_assets) != 15:
        fail(f"Expected 15 frozen V6 submission assets, found {len(publication_assets)}")
    if actual_submission_paths != expected_submission_paths:
        missing = sorted(expected_submission_paths - actual_submission_paths)
        extra = sorted(actual_submission_paths - expected_submission_paths)
        fail(f"V6 submission-asset path mismatch: missing={missing}, extra={extra}")
    role_counts = {
        role: sum(row.get("package_role") == role for row in publication_assets)
        for role in ("main_figure", "additional_file")
    }
    if role_counts != {"main_figure": 9, "additional_file": 6}:
        fail(f"Unexpected V6 submission-asset roles: {role_counts}")
    if any(len(row.get("sha256", "")) != 64 for row in publication_assets):
        fail("Invalid SHA-256 entry in PUBLICATION_ASSET_CHECKSUMS.tsv")

    with (ROOT / "CURRENT_MANUSCRIPT_ASSET_CHECKSUMS.tsv").open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        current_assets = list(csv.DictReader(handle, delimiter="\t"))
    expected_current_names = {
        "Manuscript.docx",
        *(f"Figure_{index}.png" for index in range(1, 9)),
        "Table_1.docx",
        "Table_2.docx",
        "Graphical_Abstract.png",
        "Cover_Letter.docx",
        "Additional_file_1_Supplementary_Tables_S1-S6.xlsx",
        "Additional_file_2_Supplementary_Figures_S1-S16.pdf",
        "Additional_file_3_Supplementary_Table_S7.xlsx",
        "Additional_file_4_Supplementary_Table_S8.xlsx",
        "Additional_file_5_Supplementary_Table_S10.xlsx",
        "Additional_file_6_Supplementary_Table_S9.xlsx",
    }
    actual_current_names = {row.get("filename", "") for row in current_assets}
    if len(current_assets) != 19 or actual_current_names != expected_current_names:
        missing = sorted(expected_current_names - actual_current_names)
        extra = sorted(actual_current_names - expected_current_names)
        fail(f"Current submission-asset mismatch: missing={missing}, extra={extra}")
    expected_current_roles = {
        "manuscript": 1,
        "main_figure": 8,
        "main_table": 2,
        "graphical_abstract": 1,
        "cover_letter": 1,
        "additional_file": 6,
    }
    observed_current_roles = {
        role: sum(row.get("package_role") == role for row in current_assets)
        for role in expected_current_roles
    }
    if observed_current_roles != expected_current_roles:
        fail(f"Unexpected current submission-asset roles: {observed_current_roles}")
    if any(
        len(row.get("sha256", "")) != 64 or int(row.get("bytes", "0")) <= 0
        for row in current_assets
    ):
        fail("Invalid current submission-asset size or SHA-256 entry")
    current_asset_by_name = {row["filename"]: row for row in current_assets}
    for index in range(1, 9):
        filename = f"Figure_{index}.png"
        reference = ROOT / "reference_outputs" / "current_manuscript" / filename
        row = current_asset_by_name[filename]
        data = reference.read_bytes()
        if len(data) != int(row["bytes"]) or digest_bytes(data) != row["sha256"]:
            fail(f"Current figure checksum table mismatch: {filename}")
    graphical = ROOT / "reference_outputs" / "current_manuscript" / "Graphical_Abstract.png"
    graphical_row = current_asset_by_name["Graphical_Abstract.png"]
    graphical_data = graphical.read_bytes()
    if (
        len(graphical_data) != int(graphical_row["bytes"])
        or digest_bytes(graphical_data) != graphical_row["sha256"]
    ):
        fail("Current graphical-abstract checksum table mismatch")

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
