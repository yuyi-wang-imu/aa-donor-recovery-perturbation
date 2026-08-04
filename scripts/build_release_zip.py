#!/usr/bin/env python3
"""Build and verify a byte-stable clean repository ZIP."""

from __future__ import annotations

import argparse
from pathlib import Path
import zipfile

from build_manifest import (
    OUTPUT as MANIFEST_PATH,
    ROOT,
    canonical_bytes,
    iter_release_files,
    manifest_rows,
    render_manifest,
)
from validate_repository import verify_archive


def archive_info(relative_path: str) -> zipfile.ZipInfo:
    info = zipfile.ZipInfo(relative_path, date_time=(1980, 1, 1, 0, 0, 0))
    info.create_system = 3
    info.external_attr = 0o100644 << 16
    info.compress_type = zipfile.ZIP_DEFLATED
    return info


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_zip", type=Path)
    args = parser.parse_args()
    destination = args.output_zip.resolve()
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite existing ZIP: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)

    files = list(iter_release_files())
    rows = manifest_rows()
    manifest_data = render_manifest(rows)
    if canonical_bytes(MANIFEST_PATH) != manifest_data:
        raise RuntimeError(
            "MANIFEST.tsv is stale; run scripts/build_manifest.py and "
            "scripts/validate_repository.py before packaging"
        )

    payloads = {
        path.relative_to(ROOT).as_posix(): canonical_bytes(path) for path in files
    }
    payloads["MANIFEST.tsv"] = manifest_data
    with zipfile.ZipFile(
        destination, mode="x", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as archive:
        for relative_path in sorted(payloads):
            archive.writestr(
                archive_info(relative_path),
                payloads[relative_path],
                compress_type=zipfile.ZIP_DEFLATED,
                compresslevel=9,
            )

    verify_archive(destination)
    print(f"PASS: wrote and verified {len(payloads)} files in {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
