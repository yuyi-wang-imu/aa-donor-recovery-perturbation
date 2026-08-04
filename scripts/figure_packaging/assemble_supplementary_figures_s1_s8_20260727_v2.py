"""Assemble the final eight-page supplementary-figure package.

This is a packaging-only workflow. It renumbers accepted one-page PDF assets
after excluding panels that duplicate the main figures. It does not read or
recompute analytical source data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from pypdf import PdfReader, PdfWriter


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing the eight accepted one-page PDF assets.",
    )
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    source_names = [
        "Supplementary_Figure_S1_20260726.pdf",
        "Supplementary_Figure_S2_20260726.pdf",
        "Supplementary_Figure_S3_DetectionCoverageOnly_20260726.pdf",
        "Supplementary_Figure_S4_20260726.pdf",
        "Supplementary_Figure_S5_20260726.pdf",
        "Supplementary_Figure_S6_ModuleLocalizationRobustness_20260726.pdf",
        "Supplementary_Figure_S7_DockingMatrixOnly_20260726.pdf",
        "Supplementary_Figure_S8_MatchedControlSensitivity_20260726.pdf",
    ]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    combined_path = args.output_dir / "Supplementary_Figures_S1-S8.pdf"
    archive_path = args.output_dir / "Supplementary_Figures_S1-S8_Independent_PDFs.zip"
    manifest_path = args.output_dir / "Supplementary_Figures_S1-S8_manifest.json"
    for output in (combined_path, archive_path, manifest_path):
        if output.exists():
            raise FileExistsError(f"Refusing to overwrite: {output}")

    combined = PdfWriter()
    independent_paths: list[Path] = []
    manifest: list[dict[str, str]] = []

    for number, name in enumerate(source_names, start=1):
        source = args.input_dir / name
        if not source.exists():
            raise FileNotFoundError(source)
        reader = PdfReader(str(source))
        if len(reader.pages) != 1:
            raise ValueError(f"Expected one page: {source}")

        combined.add_page(reader.pages[0])
        independent = args.output_dir / f"Supplementary_Figure_S{number}.pdf"
        writer = PdfWriter()
        writer.add_page(reader.pages[0])
        with independent.open("wb") as handle:
            writer.write(handle)
        independent_paths.append(independent)
        manifest.append(
            {
                "figure": f"Figure S{number}",
                "source_file": source.name,
                "source_sha256": sha256(source),
                "output_file": independent.name,
                "output_sha256": sha256(independent),
            }
        )

    with combined_path.open("wb") as handle:
        combined.write(handle)
    if len(PdfReader(str(combined_path)).pages) != 8:
        raise RuntimeError("Combined supplementary PDF must contain eight pages.")

    with ZipFile(archive_path, "w", ZIP_DEFLATED) as archive:
        for path in independent_paths:
            archive.write(path, arcname=path.name)

    payload = {
        "combined_pdf": combined_path.name,
        "combined_sha256": sha256(combined_path),
        "independent_archive": archive_path.name,
        "independent_archive_sha256": sha256(archive_path),
        "mapping": manifest,
    }
    manifest_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"Created {combined_path}")
    print(f"Created {archive_path}")
    print(f"Created {manifest_path}")


if __name__ == "__main__":
    main()
