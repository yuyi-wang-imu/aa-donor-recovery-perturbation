"""Versioned recovery wrapper for the v1 renderer.

The source Table C1 contains three donor-fraction rows followed by two
correlation-summary rows.  Version 1 treated all five rows as fractions and
stopped before completing the supplementary panels.  This wrapper preserves
the v1 script and changes only that input-row selection.  It also directs all
new outputs to v2 filenames so no v1 artifact is overwritten.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
V1 = HERE / "build_publication_geneformer_panels_20260803_v1.py"

spec = importlib.util.spec_from_file_location("publication_panels_v1", V1)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Could not load {V1}")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

valid_baselines = {
    "expression_centroid_shift",
    "ridge_shift",
    "geneformer_observed_recovery_shift",
}
module.baselines = module.baselines.loc[module.baselines["baseline"].isin(valid_baselines)].copy()

original_save_panel = module.save_panel
original_make_layout_proof = module.make_layout_proof


def save_panel_v2(fig, stem: str):
    return original_save_panel(fig, stem.replace("_v1", "_v2"))


def make_layout_proof_v2(stems, output_stem, ncols=2, bg="#F3F5F7"):
    stems_v2 = [stem.replace("_v1", "_v2") for stem in stems]
    return original_make_layout_proof(
        stems_v2,
        output_stem.replace("_v1", "_v2"),
        ncols=ncols,
        bg=bg,
    )


def write_source_map_v2(outputs):
    rows = []
    for panel, input_keys in outputs.items():
        for key in input_keys:
            rows.append(
                {
                    "panel": panel,
                    "source_role": key,
                    "source_path": str(module.FILES[key]),
                }
            )
    pd.DataFrame(rows).to_csv(
        module.OUT / "New_Figure_Panel_Source_Map_20260803_v2.tsv",
        sep="\t",
        index=False,
    )


module.save_panel = save_panel_v2
module.make_layout_proof = make_layout_proof_v2
module.write_source_map = write_source_map_v2


if __name__ == "__main__":
    module.main()
