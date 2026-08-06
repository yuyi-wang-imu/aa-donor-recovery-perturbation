#!/usr/bin/env python3
"""Compare a publication replay with the frozen BMC figure references.

Exact SHA-256 equality is required where the renderer is byte-stable. Figures
3 and 7 are compared in pixel space because Matplotlib font rasterization and
PNG metadata vary across supported versions even when the data and geometry
are unchanged.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image


EXACT = [1, 2, 4, 5, 6, 8, 9]
PIXEL_THRESHOLDS = {
    3: {"mean_absolute_error": 1.0, "changed_pixel_fraction": 0.20},
    7: {"mean_absolute_error": 2.5, "changed_pixel_fraction": 0.02},
}
SUPPLEMENTARY_PIXEL_THRESHOLDS = {
    "S8": {"mean_absolute_error": 0.2, "changed_pixel_fraction": 0.005},
}
SUPPLEMENTARY_EXACT = ["S9", "S10"]
SUPPLEMENTARY_PDF_SHA256 = (
    "12907c5947c5c44e32618c8d22f1ebfe5af6a9604d564310ee402a5d16335888"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pixel_metrics(actual: Path, reference: Path) -> dict[str, object]:
    with Image.open(actual) as image:
        actual_array = np.asarray(image.convert("RGB"), dtype=np.int16)
    with Image.open(reference) as image:
        reference_array = np.asarray(image.convert("RGB"), dtype=np.int16)
    if actual_array.shape != reference_array.shape:
        return {
            "pass": False,
            "reason": "dimension mismatch",
            "actual_shape": list(actual_array.shape),
            "reference_shape": list(reference_array.shape),
        }
    difference = np.abs(actual_array - reference_array)
    return {
        "mean_absolute_error": float(difference.mean()),
        "changed_pixel_fraction": float(np.any(difference, axis=2).mean()),
        "maximum_channel_error": int(difference.max()),
        "actual_shape": list(actual_array.shape),
        "reference_shape": list(reference_array.shape),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--actual-root", required=True, type=Path)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    report: dict[str, object] = {
        "policy": {
            "exact_main_figures": EXACT,
            "pixel_tolerances": PIXEL_THRESHOLDS,
            "supplementary_pixel_tolerances": SUPPLEMENTARY_PIXEL_THRESHOLDS,
            "exact_supplementary_figures": SUPPLEMENTARY_EXACT,
        },
        "results": {},
    }
    failures: list[str] = []

    for number in EXACT:
        label = f"Figure_{number}"
        actual = args.actual_root / "main_figures" / f"{label}.png"
        reference = args.reference_root / "main_figures" / f"{label}.png"
        if not actual.is_file() or not reference.is_file():
            failures.append(f"{label}: missing file")
            report["results"][label] = {"pass": False, "reason": "missing file"}
            continue
        actual_hash, reference_hash = sha256(actual), sha256(reference)
        passed = actual_hash == reference_hash
        report["results"][label] = {
            "pass": passed,
            "comparison": "exact_sha256",
            "actual_sha256": actual_hash,
            "reference_sha256": reference_hash,
        }
        if not passed:
            failures.append(f"{label}: SHA-256 mismatch")

    for number, limits in PIXEL_THRESHOLDS.items():
        label = f"Figure_{number}"
        actual = args.actual_root / "main_figures" / f"{label}.png"
        reference = args.reference_root / "main_figures" / f"{label}.png"
        if not actual.is_file() or not reference.is_file():
            failures.append(f"{label}: missing file")
            report["results"][label] = {"pass": False, "reason": "missing file"}
            continue
        metrics = pixel_metrics(actual, reference)
        passed = (
            "mean_absolute_error" in metrics
            and metrics["mean_absolute_error"] <= limits["mean_absolute_error"]
            and metrics["changed_pixel_fraction"]
            <= limits["changed_pixel_fraction"]
        )
        metrics.update(
            {
                "pass": passed,
                "comparison": "pixel_tolerance",
                "thresholds": limits,
                "actual_sha256": sha256(actual),
                "reference_sha256": sha256(reference),
            }
        )
        report["results"][label] = metrics
        if not passed:
            failures.append(f"{label}: pixel tolerance exceeded")

    for suffix in SUPPLEMENTARY_EXACT:
        label = f"Figure_{suffix}"
        actual = args.actual_root / "supplementary_figures" / f"{label}.png"
        reference = args.reference_root / "supplementary_figures" / f"{label}.png"
        if not actual.is_file() or not reference.is_file():
            failures.append(f"{label}: missing file")
            report["results"][label] = {"pass": False, "reason": "missing file"}
            continue
        actual_hash, reference_hash = sha256(actual), sha256(reference)
        passed = actual_hash == reference_hash
        report["results"][label] = {
            "pass": passed,
            "comparison": "exact_sha256",
            "actual_sha256": actual_hash,
            "reference_sha256": reference_hash,
        }
        if not passed:
            failures.append(f"{label}: SHA-256 mismatch")

    for suffix, limits in SUPPLEMENTARY_PIXEL_THRESHOLDS.items():
        label = f"Figure_{suffix}"
        actual = args.actual_root / "supplementary_figures" / f"{label}.png"
        reference = args.reference_root / "supplementary_figures" / f"{label}.png"
        if not actual.is_file() or not reference.is_file():
            failures.append(f"{label}: missing file")
            report["results"][label] = {"pass": False, "reason": "missing file"}
            continue
        metrics = pixel_metrics(actual, reference)
        passed = (
            "mean_absolute_error" in metrics
            and metrics["mean_absolute_error"] <= limits["mean_absolute_error"]
            and metrics["changed_pixel_fraction"]
            <= limits["changed_pixel_fraction"]
        )
        metrics.update(
            {
                "pass": passed,
                "comparison": "pixel_tolerance",
                "thresholds": limits,
                "actual_sha256": sha256(actual),
                "reference_sha256": sha256(reference),
            }
        )
        report["results"][label] = metrics
        if not passed:
            failures.append(f"{label}: pixel tolerance exceeded")

    supplementary_pdf = (
        args.actual_root / "supplementary_figures" / "Supplementary_Figures_S1-S8.pdf"
    )
    label = "Supplementary_Figures_S1-S8.pdf"
    if not supplementary_pdf.is_file():
        failures.append(f"{label}: missing file")
        report["results"][label] = {"pass": False, "reason": "missing file"}
    else:
        actual_hash = sha256(supplementary_pdf)
        passed = actual_hash == SUPPLEMENTARY_PDF_SHA256
        report["results"][label] = {
            "pass": passed,
            "comparison": "exact_sha256",
            "actual_sha256": actual_hash,
            "reference_sha256": SUPPLEMENTARY_PDF_SHA256,
        }
        if not passed:
            failures.append(f"{label}: SHA-256 mismatch")

    report["pass"] = not failures
    report["failures"] = failures
    payload = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
