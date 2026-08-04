"""Second versioned layout correction.

Version 3 retains every data value from v2 and only fixes two visual-QC items:
the baseline-correlation note in Figure S9A and clustered labels in Figure S10D.
All outputs are written with v3 filenames.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
V2 = HERE / "build_bmc_geneformer_panels_20260803_v2.py"

spec = importlib.util.spec_from_file_location("bmc_panels_v2", V2)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Could not load {V2}")
wrapper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wrapper)
m = wrapper.module


def save_panel_v3(fig, stem: str):
    stem = stem.replace("_v1", "_v3").replace("_v2", "_v3")
    return wrapper.original_save_panel(fig, stem)


def make_layout_proof_v3(stems, output_stem, ncols=2, bg="#F3F5F7"):
    stems = [stem.replace("_v1", "_v3").replace("_v2", "_v3") for stem in stems]
    output_stem = output_stem.replace("_v1", "_v3").replace("_v2", "_v3")
    return wrapper.original_make_layout_proof(stems, output_stem, ncols=ncols, bg=bg)


def write_source_map_v3(outputs):
    rows = []
    for panel, input_keys in outputs.items():
        for key in input_keys:
            rows.append(
                {
                    "panel": panel,
                    "source_role": key,
                    "source_path": str(m.FILES[key]),
                }
            )
    pd.DataFrame(rows).to_csv(
        m.OUT / "New_Figure_Panel_Source_Map_20260803_v3.tsv",
        sep="\t",
        index=False,
    )


def panel_s9a_v3():
    df = m.baselines.copy()
    labels = {
        "expression_centroid_shift": "Expression centroid",
        "ridge_shift": "Ridge recovery axis",
        "geneformer_observed_recovery_shift": "Geneformer embedding",
    }
    df["label"] = df["baseline"].map(labels)
    fig, ax = m.make_panel("A", "Recovery-axis baseline comparison", left=0.27, top=0.78)
    y = list(range(len(df)))[::-1]
    vals = df["positive_fraction"].to_numpy()
    ax.barh(y, vals, height=0.48, color=[m.GRAY, m.BLUE, m.PURPLE], edgecolor="none")
    for yi, val, n in zip(y, vals, df["n_donors"]):
        ax.text(min(val + 0.018, 1.04), yi, f"{int(round(val*n))}/{int(n)}", va="center", fontsize=5.5)
    ax.set_yticks(y, df["label"])
    ax.set_xlim(0, 1.08)
    ax.set_xlabel("Fraction of donors shifting toward healthy state")
    fig.text(0.27, 0.835, "Geneformer vs centroid/ridge rank correlations: 0.650 / 0.675",
             ha="left", va="bottom", fontsize=5.2, color=m.GRAY)
    m.clean_axes(ax)
    return m.save_panel(fig, "Figure_S9A_Recovery_Axis_Baseline_Comparison_20260803_v1")


def panel_s10d_v3():
    df = m.cross_model.copy()
    fig, ax = m.make_panel("D", "Complementary perturbation outputs", left=0.19)
    ax.scatter(df["geneformer_signed_deletion_shift"], df["sctenifold_significant_response_genes"],
               s=25, c=m.emphasize_colors(df["candidate"]), edgecolor="white", linewidth=0.35)
    offsets = {
        "TOP2A": (4, 3), "TERT": (4, 3), "KIT": (4, 3), "GSK3B": (4, 4),
        "PARP1": (4, 8), "CDK6": (4, -7), "CD38": (4, 2), "HIF1A": (4, -10),
        "SYK": (4, 5), "CA2": (4, 3),
    }
    for _, row in df.iterrows():
        ax.annotate(row["candidate"],
                    (row["geneformer_signed_deletion_shift"], row["sctenifold_significant_response_genes"]),
                    xytext=offsets.get(row["candidate"], (4, 3)), textcoords="offset points", fontsize=4.8)
    ax.axvline(0, color="#C9CFD7", lw=0.6)
    ax.set_xlabel("Geneformer deletion shift")
    ax.set_ylabel("scTenifoldKnk response-gene count")
    ax.text(0.98, 0.84, "Spearman ρ = -0.130", transform=ax.transAxes,
            ha="right", va="top", fontsize=5.4)
    m.clean_axes(ax, None)
    return m.save_panel(fig, "Figure_S10D_Complementary_Perturbation_Outputs_20260803_v1")


m.save_panel = save_panel_v3
m.make_layout_proof = make_layout_proof_v3
m.write_source_map = write_source_map_v3
m.panel_s9a = panel_s9a_v3
m.panel_s10d = panel_s10d_v3


if __name__ == "__main__":
    m.main()
