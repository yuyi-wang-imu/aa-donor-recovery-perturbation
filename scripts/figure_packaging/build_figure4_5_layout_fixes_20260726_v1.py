#!/usr/bin/env python3
"""Apply typography-only layout fixes to manuscript Figures 4 and 5.

The script does not recalculate, recolor, rescale, or otherwise alter any
scientific values. Figure 4 removes a verified empty spacer band between the
panel headings and their data regions. Figure 5 replaces only the embedded
A/B panel letters with larger, consistently positioned labels.
"""

from __future__ import annotations

import argparse
import pathlib

from PIL import Image, ImageDraw, ImageFont


def refuse_existing(path: pathlib.Path) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")


def fix_figure4(source: pathlib.Path, output: pathlib.Path) -> None:
    refuse_existing(output)
    image = Image.open(source).convert("RGB")
    if image.size != (4500, 4918):
        raise ValueError(f"Unexpected Figure 4 dimensions: {image.size}")
    top = image.crop((0, 0, 4500, 430))
    bottom = image.crop((0, 660, 4500, 4918))
    revised = Image.new("RGB", (4500, top.height + bottom.height), "white")
    revised.paste(top, (0, 0))
    revised.paste(bottom, (0, top.height))
    revised.save(output, dpi=(600, 600), optimize=True)


def get_font(size: int) -> ImageFont.FreeTypeFont:
    candidates = [
        pathlib.Path(r"C:\Windows\Fonts\arialbd.ttf"),
        pathlib.Path(r"C:\Windows\Fonts\calibrib.ttf"),
        pathlib.Path(r"C:\Windows\Fonts\DejaVuSans-Bold.ttf"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return ImageFont.truetype(str(candidate), size=size)
    raise FileNotFoundError("No suitable bold sans-serif font was found.")


def fix_figure5(source: pathlib.Path, output: pathlib.Path) -> None:
    refuse_existing(output)
    image = Image.open(source).convert("RGB")
    if image.size != (8040, 4214):
        raise ValueError(f"Unexpected Figure 5 dimensions: {image.size}")
    draw = ImageDraw.Draw(image)
    font = get_font(62)
    # These white patches cover only the original small panel letters.
    draw.rectangle((125, 320, 225, 440), fill="white")
    draw.rectangle((3660, 320, 3765, 440), fill="white")
    draw.text((142, 322), "A", font=font, fill="#111111")
    draw.text((3678, 322), "B", font=font, fill="#111111")
    image.save(output, dpi=(600, 600), optimize=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure4", required=True, type=pathlib.Path)
    parser.add_argument("--figure5", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fix_figure4(args.figure4, args.output_dir / "Figure_4_layout_fixed_20260726_v1.png")
    fix_figure5(args.figure5, args.output_dir / "Figure_5_panel_labels_fixed_20260726_v1.png")


if __name__ == "__main__":
    main()
