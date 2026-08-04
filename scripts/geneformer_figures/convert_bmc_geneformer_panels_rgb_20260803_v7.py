"""Convert pure-white/black-title v6 PNGs to journal-ready RGB v7 files.

Scientific content and SVG vectors are copied unchanged.  Existing versions
are never overwritten.
"""

from pathlib import Path
import shutil

from PIL import Image


ROOT = Path(__file__).resolve().parents[2] / "outputs" / "geneformer_figures"


def v7_name(path: Path) -> Path:
    return path.with_name(path.name.replace("20260803_v6", "20260803_v7"))


png_sources = sorted(ROOT.glob("*20260803_v6.png"))
svg_sources = sorted(ROOT.glob("*20260803_v6.svg"))

if len(png_sources) != 19:
    raise RuntimeError(f"Expected 19 v6 PNG files, found {len(png_sources)}")
if len(svg_sources) != 16:
    raise RuntimeError(f"Expected 16 v6 SVG files, found {len(svg_sources)}")

for source in png_sources:
    destination = v7_name(source)
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite existing file: {destination}")
    with Image.open(source) as image:
        dpi = image.info.get("dpi", (600, 600))
        rgb = Image.new("RGB", image.size, "white")
        if image.mode in {"RGBA", "LA"}:
            alpha = image.getchannel("A")
            rgb.paste(image.convert("RGB"), mask=alpha)
        else:
            rgb.paste(image.convert("RGB"))
        rgb.save(destination, format="PNG", dpi=dpi, optimize=False)

for source in svg_sources:
    destination = v7_name(source)
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite existing file: {destination}")
    shutil.copy2(source, destination)

print(f"Created {len(png_sources)} RGB PNG files and {len(svg_sources)} SVG copies as v7.")
