#!/usr/bin/env python3
"""Rebuild publication Figures 1, 2, 4, 5, and 6 from approved intermediates.

The input PNGs are author-generated, publication-facing intermediates. They
contain no raw trajectories, model weights, credentials, private design table,
or licensed database export. The operations below are deterministic layout
operations only: Figure 1 panel replacement, Figure 4 spacer removal, and
panel-letter alignment for Figures 2, 5, and 6.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


PADDING = 72
TOP_CROP_Y = 3710

PANEL_JOBS = {
    2: {
        "font_size": 70,
        "old_boxes": {
            "A": (10, 33, 43, 66),
            "B": (2173, 33, 2201, 66),
            "C": (12, 2083, 41, 2116),
            "D": (2173, 2083, 2201, 2116),
        },
        "targets": {
            "A": (10, 49),
            "B": (2173, 49),
            "C": (10, 2114),
            "D": (2173, 2114),
        },
    },
    3: {
        "font_size": 70,
        "old_boxes": {
            "A": (183, 17, 249, 83),
            "B": (2556, 17, 2611, 83),
            "C": (187, 2335, 244, 2403),
            "D": (2543, 2336, 2598, 2402),
        },
        "targets": {
            "A": (185, 146),
            "B": (2549, 146),
            "C": (185, 2463),
            "D": (2549, 2463),
        },
    },
    5: {
        "font_size": 128,
        "erase_margin": 0,
        "old_boxes": {
            "A": (142, 334, 186, 379),
            "B": (3683, 334, 3720, 379),
        },
        "targets": {
            "A": (67, 420),
            "B": (3675, 420),
        },
    },
    6: {
        "font_size": 70,
        "old_boxes": {
            "E": (22, 2934, 54, 2975),
            "F": (1462, 2934, 1491, 2975),
            "G": (2901, 2933, 2940, 2976),
        },
        "targets": {
            "E": (89, 2984),
            "F": (1487, 2984),
            "G": (2889, 2984),
        },
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def expand(box: tuple[int, int, int, int], amount: int, size: tuple[int, int]):
    x0, y0, x1, y1 = box
    width, height = size
    return (
        max(0, x0 - amount),
        max(0, y0 - amount),
        min(width, x1 + amount),
        min(height, y1 + amount),
    )


def default_font() -> Path:
    candidates = []
    windows_dir = os.environ.get("WINDIR")
    if windows_dir:
        font_dir = Path(windows_dir) / "Fonts"
        candidates.extend(
            [
                font_dir / "arialbd.ttf",
                font_dir / "calibrib.ttf",
                font_dir / "DejaVuSans-Bold.ttf",
            ]
        )
    candidates.extend(
        [
            Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
        ]
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError("No supported bold sans-serif font was found")


def require_new_output_dir(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def assemble_figure1(asset_dir: Path, destination: Path) -> None:
    source = Image.open(asset_dir / "Figure1_pre_ppt_de.png").convert("RGB")
    panels_de = Image.open(asset_dir / "Figure1_panels_de.png").convert("RGB")
    if source.size != (4320, 5020):
        raise ValueError(f"Unexpected Figure 1 source dimensions: {source.size}")
    top = source.crop((0, 0, source.width, TOP_CROP_Y))
    de_height = round(panels_de.height * source.width / panels_de.width)
    de_resized = panels_de.resize(
        (source.width, de_height), Image.Resampling.LANCZOS
    )
    composite = Image.new("RGB", (source.width, top.height + de_height), "white")
    composite.paste(top, (0, 0))
    composite.paste(de_resized, (0, top.height))
    scale = min(source.width / composite.width, source.height / composite.height)
    fitted = composite.resize(
        (round(composite.width * scale), round(composite.height * scale)),
        Image.Resampling.LANCZOS,
    )
    canvas = Image.new("RGB", source.size, "white")
    canvas.paste(
        fitted,
        ((canvas.width - fitted.width) // 2, (canvas.height - fitted.height) // 2),
    )
    canvas.save(destination, dpi=(600, 600), optimize=True)


def align_panel_letters(
    source: Path, destination: Path, figure_number: int, font_path: Path
) -> None:
    job = PANEL_JOBS[figure_number]
    original = Image.open(source)
    if original.mode != "RGB":
        raise ValueError(
            f"Figure {figure_number} intermediate must be RGB, got {original.mode}"
        )
    canvas = Image.new(
        "RGB", (original.width + 2 * PADDING, original.height + 2 * PADDING), "white"
    )
    canvas.paste(original, (PADDING, PADDING))
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype(str(font_path), job["font_size"])
    for letter, old_box in job["old_boxes"].items():
        shifted = tuple(value + PADDING for value in old_box)
        draw.rectangle(
            expand(shifted, job.get("erase_margin", 7), canvas.size), fill="white"
        )
        x, glyph_bottom = job["targets"][letter]
        bbox = draw.textbbox((0, 0), letter, font=font, anchor="lt")
        draw_y = glyph_bottom + PADDING - bbox[3]
        draw.text(
            (x + PADDING, draw_y), letter, font=font, fill="black", anchor="lt"
        )
    canvas.save(
        destination,
        format="PNG",
        dpi=original.info.get("dpi", (600, 600)),
        optimize=True,
    )


def fix_figure4(source: Path, destination: Path) -> None:
    image = Image.open(source).convert("RGB")
    if image.size != (4500, 4918):
        raise ValueError(f"Unexpected Figure 4 dimensions: {image.size}")
    top = image.crop((0, 0, 4500, 430))
    bottom = image.crop((0, 660, 4500, 4918))
    revised = Image.new("RGB", (4500, top.height + bottom.height), "white")
    revised.paste(top, (0, 0))
    revised.paste(bottom, (0, top.height))
    revised.save(destination, dpi=(600, 600), optimize=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--asset-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--font", type=Path)
    args = parser.parse_args()
    require_new_output_dir(args.output_dir)
    font_path = args.font or default_font()
    if not font_path.is_file():
        raise FileNotFoundError(font_path)

    outputs = {
        1: args.output_dir / "Figure_1.png",
        2: args.output_dir / "Figure_2.png",
        4: args.output_dir / "Figure_4.png",
        5: args.output_dir / "Figure_5.png",
        6: args.output_dir / "Figure_6.png",
    }
    assemble_figure1(args.asset_dir, outputs[1])
    align_panel_letters(
        args.asset_dir / "Figure2_pre_alignment.png", outputs[2], 2, font_path
    )
    fix_figure4(args.asset_dir / "Figure4_pre_spacer_fix.png", outputs[4])
    align_panel_letters(
        args.asset_dir / "Figure5_pre_alignment.png", outputs[5], 5, font_path
    )
    align_panel_letters(
        args.asset_dir / "Figure6_pre_alignment.png", outputs[6], 6, font_path
    )

    manifest = {
        "scope": "deterministic publication-layout reconstruction",
        "font": str(font_path),
        "outputs": {
            str(number): {
                "path": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for number, path in outputs.items()
        },
    }
    manifest_path = args.output_dir / "static_publication_figures_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
