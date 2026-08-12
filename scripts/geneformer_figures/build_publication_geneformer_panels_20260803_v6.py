"""Pure-white, black-title journal layout revision.

Version 6 preserves every scientific value, axis, annotation, and panel geometry
from the visually reviewed v4 renderer.  It changes only two presentation
properties requested for journal review: all panel letters/titles are pure
black and all composite-layout canvases are pure white.  Earlier files remain
untouched and all outputs use v6 filenames.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


HERE = Path(__file__).resolve().parent
V4 = HERE / "build_publication_geneformer_panels_20260803_v4.py"

spec = importlib.util.spec_from_file_location("publication_panels_v4", V4)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Could not load {V4}")
wrapper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wrapper)
m = wrapper.m


def make_panel_v6(letter: str, title: str, left: float = 0.19, bottom: float = 0.18,
                  right: float = 0.965, top: float = 0.82):
    fig = m.plt.figure(figsize=m.PANEL_SIZE, facecolor="white")
    ax = fig.add_axes([left, bottom, right - left, top - bottom], facecolor="white")
    fig.text(m.LETTER_X, m.TITLE_Y, letter, ha="left", va="top", fontsize=8,
             fontweight="bold", color="#000000")
    fig.text(m.TITLE_X, m.TITLE_Y, title, ha="left", va="top", fontsize=7,
             fontweight="bold", color="#000000")
    return fig, ax


def save_panel_v6(fig, stem: str):
    stem = stem.replace("_v1", "_v6").replace("_v2", "_v6").replace("_v3", "_v6").replace("_v4", "_v6")
    png = m.OUT / f"{stem}.png"
    svg = m.OUT / f"{stem}.svg"
    if png.exists() or svg.exists():
        raise FileExistsError(f"Refusing to overwrite existing v6 output: {stem}")
    fig.savefig(png, dpi=600, facecolor="white", edgecolor="white", bbox_inches=None)
    fig.savefig(svg, facecolor="white", edgecolor="white", bbox_inches=None)
    m.plt.close(fig)
    return png, svg


def make_layout_proof_v6(stems, output_stem, ncols=2, bg="#FFFFFF"):
    stems = [stem.replace("_v1", "_v6").replace("_v2", "_v6").replace("_v3", "_v6").replace("_v4", "_v6") for stem in stems]
    output_stem = output_stem.replace("_v1", "_v6").replace("_v2", "_v6").replace("_v3", "_v6").replace("_v4", "_v6")
    output = m.OUT / f"{output_stem}.png"
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing v6 output: {output}")
    imgs = [Image.open(m.OUT / f"{stem}.png").convert("RGB") for stem in stems]
    target_w = 1500
    resized = []
    for img in imgs:
        height = int(round(img.height * target_w / img.width))
        resized.append(img.resize((target_w, height), Image.Resampling.LANCZOS))
    gutter = 70
    outer = 90
    nrows = int(np.ceil(len(resized) / ncols))
    cell_h = max(img.height for img in resized)
    canvas_w = outer * 2 + ncols * target_w + (ncols - 1) * gutter
    canvas_h = outer * 2 + nrows * cell_h + (nrows - 1) * gutter
    canvas = Image.new("RGB", (canvas_w, canvas_h), bg)
    for i, img in enumerate(resized):
        row, col = divmod(i, ncols)
        x = outer + col * (target_w + gutter)
        y = outer + row * (cell_h + gutter)
        canvas.paste(img, (x, y))
    canvas.save(output, dpi=(300, 300))


def write_source_map_v6(outputs):
    rows = []
    for panel, input_keys in outputs.items():
        for key in input_keys:
            rows.append({"panel": panel, "source_role": key, "source_path": str(m.FILES[key])})
    destination = m.OUT / "New_Figure_Panel_Source_Map_20260803_v6.tsv"
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite existing v6 output: {destination}")
    pd.DataFrame(rows).to_csv(destination, sep="\t", index=False)


m.make_panel = make_panel_v6
m.save_panel = save_panel_v6
m.make_layout_proof = make_layout_proof_v6
m.write_source_map = write_source_map_v6


if __name__ == "__main__":
    m.main()
