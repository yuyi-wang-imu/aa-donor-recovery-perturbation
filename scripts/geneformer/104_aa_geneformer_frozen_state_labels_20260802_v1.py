from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


HSPC = ["CD34", "KIT", "PROM1", "AVP", "GATA2", "MECOM", "MEIS1", "HLF", "SPINK2"]
MEGAKARYOCYTE = ["PPBP", "PF4", "ITGA2B", "GP9", "MPL", "VWF"]
MIN_BEST_SCORE = 0.02
MIN_MARGIN = 0.005


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def score_set(log_cp10k: sp.csr_matrix, genes: list[str], gene_to_idx: dict[str, int]) -> np.ndarray:
    idx = [gene_to_idx[g] for g in genes if g in gene_to_idx]
    if not idx:
        return np.full(log_cp10k.shape[0], np.nan)
    return np.asarray(log_cp10k[:, idx].mean(axis=1)).ravel()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    input_dir, output_dir = Path(args.input_dir), Path(args.output_dir)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise RuntimeError(f"Refusing non-empty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = input_dir / "GSE247531_CD34_balanced_cells_by_genes_mvp_v1.mtx"
    genes_path = input_dir / "GSE247531_CD34_balanced_gene_symbols_mvp_v1.tsv"
    meta_path = input_dir / "GSE247531_CD34_balanced_cell_metadata_mvp_v1.tsv"
    x = scipy.io.mmread(matrix_path).tocsr().astype(np.float64)
    genes = pd.read_csv(genes_path, sep="\t")["gene_symbol"].astype(str).tolist()
    meta = pd.read_csv(meta_path, sep="\t")
    if x.shape != (len(meta), len(genes)):
        raise RuntimeError(f"Matrix mismatch: {x.shape}, meta={len(meta)}, genes={len(genes)}")
    lib = np.asarray(x.sum(axis=1)).ravel()
    if np.any(lib <= 0):
        raise RuntimeError("Zero-library cells are not allowed")
    log_cp10k = x.multiply((10000.0 / lib)[:, None]).tocsr()
    log_cp10k.data = np.log1p(log_cp10k.data)
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    hspc = score_set(log_cp10k, HSPC, gene_to_idx)
    mega = score_set(log_cp10k, MEGAKARYOCYTE, gene_to_idx)
    best = np.maximum(hspc, mega)
    margin = np.abs(hspc - mega)
    labels = np.where(hspc >= mega, "HSPC-marker-class", "megakaryocyte-marker-class")
    labels[(best < MIN_BEST_SCORE) | (margin < MIN_MARGIN)] = "Unclassified"
    out = meta.copy()
    out["hspc_marker_score"] = hspc
    out["megakaryocyte_marker_score"] = mega
    out["marker_best_score"] = best
    out["marker_margin"] = margin
    out["frozen_state_label"] = labels
    out.to_csv(output_dir / "frozen_state_labels.tsv", sep="\t", index=False)
    counts = out.groupby(["subject", "disease", "timepoint", "frozen_state_label"], dropna=False).size().rename("n_cells").reset_index()
    counts.to_csv(output_dir / "state_counts_by_subject_timepoint.tsv", sep="\t", index=False)
    summary = out["frozen_state_label"].value_counts(dropna=False).rename_axis("frozen_state_label").reset_index(name="n_cells")
    summary["fraction"] = summary.n_cells / len(out)
    summary.to_csv(output_dir / "state_label_summary.tsv", sep="\t", index=False)
    qc = {
        "matrix_shape": list(x.shape), "n_cells": int(len(out)), "n_genes": int(len(genes)),
        "marker_sets": {"HSPC": HSPC, "Megakaryocyte": MEGAKARYOCYTE},
        "present_markers": {"HSPC": [g for g in HSPC if g in gene_to_idx], "Megakaryocyte": [g for g in MEGAKARYOCYTE if g in gene_to_idx]},
        "min_best_score": MIN_BEST_SCORE, "min_margin": MIN_MARGIN,
        "label_counts": {str(k): int(v) for k, v in out.frozen_state_label.value_counts().items()},
        "all_scores_finite": bool(np.isfinite(hspc).all() and np.isfinite(mega).all()),
        "input_sha256": {str(p): sha256(p) for p in [matrix_path, genes_path, meta_path]},
        "classification_note": "Predefined marker scoring only; no clustering, perturbation-result tuning, or threshold optimization.",
    }
    (output_dir / "technical_qc.json").write_text(json.dumps(qc, indent=2, ensure_ascii=False), encoding="utf-8")
    print("STATE_LABELS_COMPLETE")
    print(json.dumps(qc["label_counts"], ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
