"""Final visual-QC candidate renderer.

Version 4 retains all values and geometry from v3 and separates the dense
low-response labels in Figure S10D using data-coordinate offsets.  Outputs are
written with v4 filenames; earlier versions remain untouched.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
V3 = HERE / "build_bmc_geneformer_panels_20260803_v3.py"

spec = importlib.util.spec_from_file_location("bmc_panels_v3", V3)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Could not load {V3}")
wrapper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wrapper)
m = wrapper.m


def save_panel_v4(fig, stem: str):
    stem = stem.replace("_v1", "_v4").replace("_v2", "_v4").replace("_v3", "_v4")
    return wrapper.wrapper.original_save_panel(fig, stem)


def make_layout_proof_v4(stems, output_stem, ncols=2, bg="#F3F5F7"):
    stems = [stem.replace("_v1", "_v4").replace("_v2", "_v4").replace("_v3", "_v4") for stem in stems]
    output_stem = output_stem.replace("_v1", "_v4").replace("_v2", "_v4").replace("_v3", "_v4")
    return wrapper.wrapper.original_make_layout_proof(stems, output_stem, ncols=ncols, bg=bg)


def write_source_map_v4(outputs):
    rows = []
    for panel, input_keys in outputs.items():
        for key in input_keys:
            rows.append({"panel": panel, "source_role": key, "source_path": str(m.FILES[key])})
    pd.DataFrame(rows).to_csv(
        m.OUT / "New_Figure_Panel_Source_Map_20260803_v4.tsv",
        sep="\t",
        index=False,
    )


def panel_s10d_v4():
    df = m.cross_model.copy()
    fig, ax = m.make_panel("D", "Complementary perturbation outputs", left=0.19)
    ax.scatter(df["geneformer_signed_deletion_shift"], df["sctenifold_significant_response_genes"],
               s=25, c=m.emphasize_colors(df["candidate"]), edgecolor="white", linewidth=0.35)
    y_offsets = {
        "TOP2A": 2.3, "TERT": 2.3, "KIT": 2.3, "GSK3B": 2.3,
        "PARP1": 6.0, "CDK6": -2.0, "CD38": 2.0, "HIF1A": -3.5,
        "SYK": 4.0, "CA2": 2.0,
    }
    for _, row in df.iterrows():
        x = float(row["geneformer_signed_deletion_shift"])
        y = float(row["sctenifold_significant_response_genes"])
        ax.text(x + 0.0012, y + y_offsets.get(row["candidate"], 2.0), row["candidate"],
                fontsize=4.8, ha="left", va="center")
    ax.axvline(0, color="#C9CFD7", lw=0.6)
    ax.set_ylim(4, 105)
    ax.set_xlabel("Geneformer deletion shift")
    ax.set_ylabel("scTenifoldKnk response-gene count")
    ax.text(0.98, 0.84, "Spearman ρ = -0.130", transform=ax.transAxes,
            ha="right", va="top", fontsize=5.4)
    m.clean_axes(ax, None)
    return m.save_panel(fig, "Figure_S10D_Complementary_Perturbation_Outputs_20260803_v1")


m.save_panel = save_panel_v4
m.make_layout_proof = make_layout_proof_v4
m.write_source_map = write_source_map_v4
m.panel_s10d = panel_s10d_v4


if __name__ == "__main__":
    m.main()
