#!/usr/bin/env python3
"""Create a path-neutral SHA-256 manifest for canonical release bytes."""

from __future__ import annotations

import csv
import hashlib
import io
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "MANIFEST.tsv"
EXCLUDED_DIRS = {"__pycache__", ".git", ".Rproj.user", "outputs"}
EXCLUDED_NAMES = {
    "MANIFEST.tsv",
    ".DS_Store",
    "Thumbs.db",
    ".RData",
    ".Rhistory",
    "UPLOAD_READINESS_AUDIT_20260804_v1.md",
    "PROPOSED_UPLOAD_MANIFEST_20260804_v1.tsv",
}
EXCLUDED_SUFFIXES = {".pyc", ".pyo"}


def canonical_bytes(path: Path) -> bytes:
    """Return the bytes stored in release archives.

    Git's text=auto/eol=lf policy normalizes CRLF to LF for UTF-8 text while
    leaving binary content unchanged. Applying the same rule here keeps
    manifests reproducible even in an older Windows checkout whose worktree
    still contains CRLF files.
    """

    data = path.read_bytes()
    try:
        data.decode("utf-8")
    except UnicodeDecodeError:
        return data
    return data.replace(b"\r\n", b"\n")


def digest_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest().upper()


def iter_release_files() -> Iterable[Path]:
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT)
        if (
            any(part in EXCLUDED_DIRS for part in relative.parts)
            or path.name in EXCLUDED_NAMES
            or path.suffix.lower() in EXCLUDED_SUFFIXES
        ):
            continue
        yield path


def manifest_rows() -> list[tuple[str, int, str]]:
    rows = []
    for path in iter_release_files():
        data = canonical_bytes(path)
        rows.append(
            (
                path.relative_to(ROOT).as_posix(),
                len(data),
                digest_bytes(data),
            )
        )
    return rows


def render_manifest(rows: Iterable[tuple[str, int, str]]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(["relative_path", "bytes", "sha256"])
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def main() -> None:
    rows = manifest_rows()
    OUTPUT.write_bytes(render_manifest(rows))
    print(f"Wrote {len(rows)} entries to {OUTPUT}")


if __name__ == "__main__":
    main()
