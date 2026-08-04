from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import pandas as pd


POSITIVE_CONTROLS = ["MPL", "JAK2", "STAT5A", "STAT5B", "PIM1", "BCL2L1"]


def load_engine(path: Path):
    spec = importlib.util.spec_from_file_location("frozen_geneformer_engine", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import engine: {path}")
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return getattr(loaded, "module", loaded)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--engine", required=True)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    engine = load_engine(Path(args.engine))
    engine.MVP_TARGETS = list(POSITIVE_CONTROLS)
    engine.CANDIDATES = list(POSITIVE_CONTROLS)

    def select_positive_genes(manifest, symbol_to_token, detected_sets, meta):
        baseline_idx = meta.index[meta.analysis_group.eq("SAA_baseline")].to_numpy()
        selected, rows = [], []
        for gene in POSITIVE_CONTROLS:
            token = symbol_to_token.get(gene)
            detected_n = sum(token in detected_sets[i] for i in baseline_idx) if token is not None else 0
            rows.append({"candidate": gene, "control_gene": gene, "frozen_order": 0,
                         "detected_baseline_cells": detected_n,
                         "status": "accepted" if token is not None else "token_missing",
                         "reason": "predefined_hematopoietic_positive_control" if token is not None else "token_missing"})
            if token is None:
                raise RuntimeError(f"Positive-control gene lacks Geneformer token: {gene}")
            selected.append((gene, gene, "candidate", 0))
        return selected, pd.DataFrame(rows)

    engine.select_mvp_genes = select_positive_genes
    sys.argv = [Path(args.engine).name, "--input-dir", args.input_dir,
                "--output-dir", args.output_dir, "--mode", "mvp"]
    return int(engine.main())


if __name__ == "__main__":
    raise SystemExit(main())
