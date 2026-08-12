from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import pickle
import platform
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import scipy.io
import scipy.sparse as sp
from scipy.stats import spearmanr
import torch
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from transformers import BertForMaskedLM


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = Path(os.environ.get("GENEFORMER_SOURCE_DIR", "__GENEFORMER_SOURCE_DIR_NOT_SET__"))
MODEL_DIR = Path(os.environ.get("GENEFORMER_MODEL_DIR", "__GENEFORMER_MODEL_DIR_NOT_SET__"))
CONTROL_ROOT = REPOSITORY_ROOT / "derived_data" / "computational_perturbation"
CANDIDATES = ["CDK6", "CA2", "PARP1", "KIT", "SYK", "GSK3B", "HIF1A", "TOP2A", "TERT", "CD38"]
MVP_TARGETS = ["KIT", "SYK"]
SEED = 20260802
BOOTSTRAP_REPS = 2000
MAX_CELLS_PER_DONOR_GENE = 16
MIN_CELLS_PER_DONOR_GENE = 2
MAX_CONTROLS_PER_TARGET = 5
PERTURBATION = "overexpress"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def refuse_nonempty(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise RuntimeError(f"Refusing to write into non-empty output directory: {path}")
    path.mkdir(parents=True, exist_ok=True)


def find_one(pattern: str) -> Path:
    # The frozen source snapshot contains a build/lib mirror.  Use the
    # authoritative source-tree asset and ignore the generated mirror.
    hits = [p for p in sorted(SOURCE_DIR.rglob(pattern)) if "/build/lib/" not in p.as_posix()]
    if len(hits) != 1:
        raise RuntimeError(f"Expected exactly one non-build source asset for {pattern}; observed {len(hits)}: {hits}")
    return hits[0]


def load_pickle(path: Path):
    with path.open("rb") as handle:
        return pickle.load(handle)


def canonical_ensembl(value):
    if isinstance(value, str):
        return value.split(".")[0]
    if isinstance(value, (list, tuple, np.ndarray)) and len(value):
        return str(value[0]).split(".")[0]
    return None


def normalize_labels(meta: pd.DataFrame) -> pd.DataFrame:
    out = meta.copy()
    disease = out["disease"].astype(str).str.strip().str.lower()
    timepoint = out["timepoint"].astype(str).str.strip().str.lower()
    out["analysis_group"] = "other"
    out.loc[disease.isin(["hd", "healthy", "healthy control", "control"]) | timepoint.eq("hd"), "analysis_group"] = "HD"
    baseline = timepoint.str.contains("base|pre|pretreat|0m|0 month", regex=True)
    month6 = timepoint.str.contains("6m|6 month|six", regex=True)
    saa = disease.str.contains("saa|aplastic") | ~out["analysis_group"].eq("HD")
    out.loc[saa & baseline, "analysis_group"] = "SAA_baseline"
    out.loc[saa & month6, "analysis_group"] = "SAA_6M"
    return out


def build_collapsed_matrix(
    x: sp.csr_matrix,
    symbols: list[str],
    name_to_id: dict,
    median_dict: dict,
    token_dict: dict,
):
    old_indices = []
    ens_ids = []
    new_index: dict[str, int] = {}
    for i, symbol in enumerate(symbols):
        raw = symbol if symbol.startswith("ENSG") else name_to_id.get(symbol)
        ens = canonical_ensembl(raw)
        if not ens or ens not in median_dict or ens not in token_dict:
            continue
        if ens not in new_index:
            new_index[ens] = len(ens_ids)
            ens_ids.append(ens)
        old_indices.append((i, new_index[ens]))
    if not old_indices:
        raise RuntimeError("No input genes mapped into the Geneformer V2 vocabulary")
    rows = np.array([i for i, _ in old_indices], dtype=np.int64)
    cols = np.array([j for _, j in old_indices], dtype=np.int64)
    selector = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(x.shape[1], len(ens_ids)))
    collapsed = (x @ selector).tocsr()
    return collapsed, ens_ids


def tokenize_matrix(x: sp.csr_matrix, ens_ids: list[str], median_dict: dict, token_dict: dict):
    medians = np.asarray([float(median_dict[g]) for g in ens_ids], dtype=np.float64)
    tokens = np.asarray([int(token_dict[g]) for g in ens_ids], dtype=np.int64)
    cls_token = int(token_dict["<cls>"])
    eos_token = int(token_dict["<eos>"])
    max_genes = 4094
    sequences: list[list[int]] = []
    detected_sets: list[set[int]] = []
    token_ranks: list[dict[int, float]] = []
    for row in range(x.shape[0]):
        start, end = x.indptr[row], x.indptr[row + 1]
        idx = x.indices[start:end]
        vals = x.data[start:end].astype(np.float64, copy=False)
        scaled = vals / medians[idx]
        order = np.argsort(-scaled, kind="mergesort")[:max_genes]
        ranked = tokens[idx[order]].tolist()
        seq = [cls_token] + ranked + [eos_token]
        sequences.append(seq)
        detected_sets.append(set(ranked))
        denom = max(1, len(ranked) - 1)
        token_ranks.append({token: rank / denom for rank, token in enumerate(ranked)})
    return sequences, detected_sets, token_ranks


def embed_sequences(model, sequences: list[list[int]], pad_token: int, device: torch.device, batch_size: int = 16):
    outputs = np.empty((len(sequences), model.config.hidden_size), dtype=np.float32)
    with torch.inference_mode():
        for start in range(0, len(sequences), batch_size):
            batch = sequences[start : start + batch_size]
            width = max(len(seq) for seq in batch)
            ids = torch.full((len(batch), width), pad_token, dtype=torch.long, device=device)
            mask = torch.zeros((len(batch), width), dtype=torch.long, device=device)
            for i, seq in enumerate(batch):
                ids[i, : len(seq)] = torch.as_tensor(seq, dtype=torch.long, device=device)
                mask[i, : len(seq)] = 1
            with torch.autocast(device_type="cuda", dtype=torch.float16, enabled=device.type == "cuda"):
                result = model.bert(input_ids=ids, attention_mask=mask, return_dict=True)
                cls = result.last_hidden_state[:, 0, :]
            outputs[start : start + len(batch)] = cls.float().cpu().numpy()
    return outputs


def overexpress_sequence_pair(
    sequence: list[int], token: int, cls_token: int, eos_token: int, max_len: int
) -> tuple[list[int], list[int], int]:
    """Return official-style perturbed and overflow-matched original V2 sequences."""
    if len(sequence) < 2 or sequence[0] != cls_token or sequence[-1] != eos_token:
        raise RuntimeError("Unexpected V2 special-token layout")
    original_middle = list(sequence[1:-1])
    perturbed_middle = [token, *[value for value in original_middle if value != token]]
    overflow = max(0, len(perturbed_middle) + 2 - max_len)
    if overflow:
        perturbed_middle = perturbed_middle[:-overflow]
        original_middle = original_middle[:-overflow]
    perturbed = [cls_token, *perturbed_middle, eos_token]
    comparison = [cls_token, *original_middle, eos_token]
    if len(perturbed) > max_len or perturbed[0] != cls_token or perturbed[-1] != eos_token:
        raise RuntimeError("Overexpression sequence contract failed")
    return perturbed, comparison, overflow


def matched_original_embeddings(
    model,
    indices: list[int],
    comparison_sequences: list[list[int]],
    original_sequences: list[list[int]],
    cached_embeddings: np.ndarray,
    pad_token: int,
    device: torch.device,
) -> tuple[np.ndarray, int]:
    """Reuse cached CLS embeddings unless official overflow matching changed input."""
    output = cached_embeddings[indices].copy()
    changed = [position for position, (index, seq) in enumerate(zip(indices, comparison_sequences)) if seq != original_sequences[index]]
    if changed:
        changed_embeddings = embed_sequences(model, [comparison_sequences[position] for position in changed], pad_token, device, batch_size=16)
        output[changed] = changed_embeddings
    return output, len(changed)


def unit(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if not np.isfinite(norm) or norm <= 0:
        raise RuntimeError("Zero or non-finite recovery-axis norm")
    return vector / norm


def donor_axes(emb: np.ndarray, meta: pd.DataFrame):
    hd_idx = np.flatnonzero(meta.analysis_group.eq("HD").to_numpy())
    baseline_idx = np.flatnonzero(meta.analysis_group.eq("SAA_baseline").to_numpy())
    if len(hd_idx) == 0 or len(baseline_idx) == 0:
        raise RuntimeError(f"Missing HD or baseline cells: HD={len(hd_idx)} baseline={len(baseline_idx)}")
    hd_centroid = emb[hd_idx].mean(axis=0)
    donors = sorted(set(meta.loc[meta.analysis_group.eq("SAA_baseline"), "subject"]) & set(meta.loc[meta.analysis_group.eq("SAA_6M"), "subject"]))
    axes = {}
    observed = []
    for donor in donors:
        other_base = np.flatnonzero((meta.analysis_group.eq("SAA_baseline") & meta.subject.ne(donor)).to_numpy())
        own_base = np.flatnonzero((meta.analysis_group.eq("SAA_baseline") & meta.subject.eq(donor)).to_numpy())
        own_6m = np.flatnonzero((meta.analysis_group.eq("SAA_6M") & meta.subject.eq(donor)).to_numpy())
        if not len(other_base) or not len(own_base) or not len(own_6m):
            continue
        axis = unit(hd_centroid - emb[other_base].mean(axis=0))
        shift = float(np.dot(emb[own_6m].mean(axis=0) - emb[own_base].mean(axis=0), axis))
        axes[donor] = axis
        observed.append({"subject": donor, "geneformer_observed_recovery_shift": shift, "n_baseline_cells": len(own_base), "n_6m_cells": len(own_6m)})
    return axes, pd.DataFrame(observed)


def select_mvp_genes(manifest: pd.DataFrame, symbol_to_token: dict[str, int], detected_sets: list[set[int]], meta: pd.DataFrame):
    baseline_idx = np.flatnonzero(meta.analysis_group.eq("SAA_baseline").to_numpy())
    rows = []
    selected = []
    for target in MVP_TARGETS:
        if target not in symbol_to_token:
            raise RuntimeError(f"MVP target lacks Geneformer token: {target}")
        selected.append((target, target, "candidate", 0))
        subset = manifest[(manifest["run_role"] == "matched_control") & (manifest["candidate"] == target)].copy()
        subset = subset.sort_values("run_id", kind="stable")
        accepted = 0
        for order, record in enumerate(subset.itertuples(index=False), start=1):
            gene = str(record.gKO)
            token = symbol_to_token.get(gene)
            detected_n = sum(token in detected_sets[i] for i in baseline_idx) if token is not None else 0
            status = "accepted" if token is not None and detected_n >= 10 and accepted < MAX_CONTROLS_PER_TARGET else "skipped"
            reason = "accepted_frozen_order" if status == "accepted" else ("token_missing" if token is None else ("insufficient_detected_cells" if detected_n < 10 else "quota_reached"))
            rows.append({"candidate": target, "control_gene": gene, "frozen_order": order, "detected_baseline_cells": detected_n, "status": status, "reason": reason})
            if status == "accepted":
                selected.append((gene, target, "matched_control", order))
                accepted += 1
            if accepted == MAX_CONTROLS_PER_TARGET:
                break
        if accepted < 4:
            raise RuntimeError(f"Fewer than four tokenizable controls for {target}: {accepted}")
    return selected, pd.DataFrame(rows)


def bootstrap_summary(values: np.ndarray, rng: np.random.Generator):
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return {"n_donors": 0, "mean": np.nan, "median": np.nan, "ci_low": np.nan, "ci_high": np.nan, "positive_fraction": np.nan}
    boot = np.empty(BOOTSTRAP_REPS, dtype=float)
    for i in range(BOOTSTRAP_REPS):
        boot[i] = np.mean(rng.choice(values, size=len(values), replace=True))
    return {
        "n_donors": len(values),
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "ci_low": float(np.quantile(boot, 0.025)),
        "ci_high": float(np.quantile(boot, 0.975)),
        "positive_fraction": float(np.mean(values > 0)),
    }


def make_pseudobulk(x: sp.csr_matrix, meta: pd.DataFrame):
    labels = meta["subject"].astype(str) + "__" + meta["analysis_group"].astype(str)
    records = []
    arrays = []
    for label in sorted(labels.unique()):
        idx = np.flatnonzero(labels.eq(label).to_numpy())
        mean_counts = np.asarray(x[idx].mean(axis=0)).ravel()
        norm = np.log1p(mean_counts / max(mean_counts.sum(), 1.0) * 1e6)
        subject, group = label.split("__", 1)
        arrays.append(norm)
        records.append({"subject": subject, "analysis_group": group, "n_cells": len(idx)})
    return np.vstack(arrays), pd.DataFrame(records)


def simple_baselines(x: sp.csr_matrix, meta: pd.DataFrame, gf_observed: pd.DataFrame):
    pb, pb_meta = make_pseudobulk(x, meta)
    donors = sorted(set(pb_meta.loc[pb_meta.analysis_group.eq("SAA_baseline"), "subject"]) & set(pb_meta.loc[pb_meta.analysis_group.eq("SAA_6M"), "subject"]))
    records = []
    for donor in donors:
        hd = np.flatnonzero(pb_meta.analysis_group.eq("HD").to_numpy())
        train_base = np.flatnonzero((pb_meta.analysis_group.eq("SAA_baseline") & pb_meta.subject.ne(donor)).to_numpy())
        own_base = np.flatnonzero((pb_meta.analysis_group.eq("SAA_baseline") & pb_meta.subject.eq(donor)).to_numpy())
        own_6m = np.flatnonzero((pb_meta.analysis_group.eq("SAA_6M") & pb_meta.subject.eq(donor)).to_numpy())
        if not len(hd) or not len(train_base) or len(own_base) != 1 or len(own_6m) != 1:
            continue
        centroid_axis = unit(pb[hd].mean(axis=0) - pb[train_base].mean(axis=0))
        centroid_shift = float(np.dot(pb[own_6m[0]] - pb[own_base[0]], centroid_axis))
        train_idx = np.concatenate([hd, train_base])
        y = np.concatenate([np.ones(len(hd)), np.zeros(len(train_base))])
        scaler = StandardScaler(with_mean=True, with_std=True)
        train_scaled = scaler.fit_transform(pb[train_idx])
        ridge = Ridge(alpha=1.0).fit(train_scaled, y)
        ridge_shift = float(ridge.predict(scaler.transform(pb[own_6m]))[0] - ridge.predict(scaler.transform(pb[own_base]))[0])
        records.append({"subject": donor, "expression_centroid_shift": centroid_shift, "ridge_shift": ridge_shift})
    table = pd.DataFrame(records).merge(gf_observed[["subject", "geneformer_observed_recovery_shift"]], on="subject", how="inner")
    summary = []
    for col in ["expression_centroid_shift", "ridge_shift", "geneformer_observed_recovery_shift"]:
        vals = table[col].to_numpy(float)
        summary.append({"baseline": col, "n_donors": len(vals), "median_shift": float(np.median(vals)), "positive_fraction": float(np.mean(vals > 0))})
    for col in ["expression_centroid_shift", "ridge_shift"]:
        rho = spearmanr(table[col], table["geneformer_observed_recovery_shift"]).statistic
        summary.append({"baseline": f"spearman_{col}_vs_geneformer", "n_donors": len(table), "median_shift": float(rho), "positive_fraction": np.nan})
    return table, pd.DataFrame(summary)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--mode", choices=["mvp", "full"], default="mvp")
    args = parser.parse_args()
    if args.mode != "mvp":
        raise RuntimeError("This frozen script implements the MVP only; full analysis requires a separately versioned script")

    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    if torch.cuda.is_available(): torch.cuda.manual_seed_all(SEED)

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    refuse_nonempty(output_dir)

    matrix_path = input_dir / "GSE247531_CD34_balanced_cells_by_genes_mvp_v1.mtx"
    genes_path = input_dir / "GSE247531_CD34_balanced_gene_symbols_mvp_v1.tsv"
    meta_path = input_dir / "GSE247531_CD34_balanced_cell_metadata_mvp_v1.tsv"
    manifest_path = CONTROL_ROOT / "run_manifest_210.csv"
    for path in [matrix_path, genes_path, meta_path, manifest_path, MODEL_DIR / "model.safetensors"]:
        if not path.exists(): raise FileNotFoundError(path)

    token_path = find_one("token_dictionary_gc104M.pkl")
    median_path = find_one("gene_median_dictionary_gc104M.pkl")
    name_id_path = find_one("gene_name_id_dict_gc104M.pkl")
    token_dict = load_pickle(token_path)
    median_dict = load_pickle(median_path)
    name_to_id = load_pickle(name_id_path)
    symbol_to_token = {symbol: token_dict.get(canonical_ensembl(ens)) for symbol, ens in name_to_id.items()}
    symbol_to_token = {k: int(v) for k, v in symbol_to_token.items() if v is not None}

    x = scipy.io.mmread(matrix_path).tocsr().astype(np.float32)
    genes = pd.read_csv(genes_path, sep="\t")["gene_symbol"].astype(str).tolist()
    meta = normalize_labels(pd.read_csv(meta_path, sep="\t"))
    if x.shape != (len(meta), len(genes)):
        raise RuntimeError(f"Matrix/metadata mismatch: {x.shape}, {len(meta)}, {len(genes)}")
    collapsed, ens_ids = build_collapsed_matrix(x, genes, name_to_id, median_dict, token_dict)
    sequences, detected_sets, token_ranks = tokenize_matrix(collapsed, ens_ids, median_dict, token_dict)

    group_counts = meta.analysis_group.value_counts().to_dict()
    if group_counts.get("HD", 0) < 64 or group_counts.get("SAA_baseline", 0) < 64 or group_counts.get("SAA_6M", 0) < 64:
        raise RuntimeError(f"Insufficient analysis groups after normalization: {group_counts}")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = BertForMaskedLM.from_pretrained(str(MODEL_DIR), local_files_only=True, torch_dtype=torch.float16 if device.type == "cuda" else torch.float32)
    model.eval().to(device)
    pad_token = int(token_dict["<pad>"])
    embeddings = embed_sequences(model, sequences, pad_token, device, batch_size=16)
    axes, observed = donor_axes(embeddings, meta)
    observed.to_csv(output_dir / "recovery_axis_lodo.tsv", sep="\t", index=False)

    run_manifest = pd.read_csv(manifest_path)
    selected, selection_audit = select_mvp_genes(run_manifest, symbol_to_token, detected_sets, meta)
    selection_audit.to_csv(output_dir / "matched_control_selection_audit.tsv", sep="\t", index=False)
    pd.DataFrame(selected, columns=["gene", "candidate", "run_role", "frozen_order"]).to_csv(output_dir / "selected_genes.tsv", sep="\t", index=False)

    donor_rows = []
    donor_rows_all_baseline = []
    metric_rows = []
    rng = np.random.default_rng(SEED)
    baseline_mask = meta.analysis_group.eq("SAA_baseline").to_numpy()
    baseline_indices = np.flatnonzero(baseline_mask)
    model_max_len = int(model.config.max_position_embeddings)
    overflow_matched_comparisons = 0
    maximum_perturbed_length = 0
    for gene, candidate, role, frozen_order in selected:
        token = symbol_to_token[gene]
        expressing = [i for i in baseline_indices if token in detected_sets[i]]
        expression_fraction = len(expressing) / max(1, len(baseline_indices))
        ranks = [token_ranks[i][token] for i in expressing]
        for donor, axis in axes.items():
            # Secondary, official Geneformer overexpression scope: all baseline
            # cells, including cells in which the target token was not detected.
            # This block must not depend on the detected-cell primary endpoint.
            all_donor_cells = [i for i in baseline_indices if str(meta.iloc[i]["subject"]) == str(donor)]
            if len(all_donor_cells) >= MIN_CELLS_PER_DONOR_GENE:
                all_rng = np.random.default_rng(SEED + 100000 + sum(ord(c) for c in gene) + sum(ord(c) for c in str(donor)))
                if len(all_donor_cells) > MAX_CELLS_PER_DONOR_GENE:
                    all_donor_cells = sorted(all_rng.choice(all_donor_cells, MAX_CELLS_PER_DONOR_GENE, replace=False).tolist())
                all_pairs = [overexpress_sequence_pair(sequences[i], token, int(token_dict["<cls>"]), int(token_dict["<eos>"]), model_max_len) for i in all_donor_cells]
                overexpressed_all = [pair[0] for pair in all_pairs]
                comparison_all = [pair[1] for pair in all_pairs]
                overflow_matched_comparisons += sum(pair[2] > 0 for pair in all_pairs)
                maximum_perturbed_length = max(maximum_perturbed_length, max(map(len, overexpressed_all)))
                overexpressed_all_emb = embed_sequences(model, overexpressed_all, pad_token, device, batch_size=16)
                original_all_emb, _ = matched_original_embeddings(model, all_donor_cells, comparison_all, sequences, embeddings, pad_token, device)
                all_shifts = (overexpressed_all_emb - original_all_emb) @ axis
                donor_rows_all_baseline.append({"gene": gene, "candidate": candidate, "run_role": role, "subject": donor, "n_cells": len(all_donor_cells), "cell_scope": "all_baseline_cells_official_overexpression", "mean_overexpression_recovery_shift": float(np.mean(all_shifts)), "median_overexpression_recovery_shift": float(np.median(all_shifts))})

            # Primary endpoint: the exact detected-cell scope used by the frozen
            # deletion analysis, enabling a paired bidirectional comparison.
            donor_cells = [i for i in expressing if str(meta.iloc[i]["subject"]) == str(donor)]
            if len(donor_cells) < MIN_CELLS_PER_DONOR_GENE:
                continue
            donor_rng = np.random.default_rng(SEED + sum(ord(c) for c in gene) + sum(ord(c) for c in str(donor)))
            if len(donor_cells) > MAX_CELLS_PER_DONOR_GENE:
                donor_cells = sorted(donor_rng.choice(donor_cells, MAX_CELLS_PER_DONOR_GENE, replace=False).tolist())
            detected_pairs = [overexpress_sequence_pair(sequences[i], token, int(token_dict["<cls>"]), int(token_dict["<eos>"]), model_max_len) for i in donor_cells]
            overexpressed = [pair[0] for pair in detected_pairs]
            comparison_detected = [pair[1] for pair in detected_pairs]
            overflow_matched_comparisons += sum(pair[2] > 0 for pair in detected_pairs)
            maximum_perturbed_length = max(maximum_perturbed_length, max(map(len, overexpressed)))
            overexpressed_emb = embed_sequences(model, overexpressed, pad_token, device, batch_size=16)
            original_emb, _ = matched_original_embeddings(model, donor_cells, comparison_detected, sequences, embeddings, pad_token, device)
            shifts = (overexpressed_emb - original_emb) @ axis
            donor_rows.append({"gene": gene, "candidate": candidate, "run_role": role, "subject": donor, "n_cells": len(donor_cells), "cell_scope": "detected_same_cells_as_deletion", "mean_overexpression_recovery_shift": float(np.mean(shifts)), "median_overexpression_recovery_shift": float(np.median(shifts))})
        metric_rows.append({"gene": gene, "candidate": candidate, "run_role": role, "frozen_order": frozen_order, "baseline_detection_fraction": expression_fraction, "median_normalized_token_rank": float(np.median(ranks)) if ranks else np.nan, "n_expressing_baseline_cells": len(expressing)})

    donor_table = pd.DataFrame(donor_rows)
    donor_table.to_csv(output_dir / "donor_gene_effects.tsv", sep="\t", index=False)
    donor_all_table = pd.DataFrame(donor_rows_all_baseline)
    donor_all_table.to_csv(output_dir / "donor_gene_effects_all_baseline.tsv", sep="\t", index=False)
    if donor_table.empty:
        raise RuntimeError("No donor-level perturbation effects were generated")
    if donor_all_table.empty:
        raise RuntimeError("No all-baseline overexpression effects were generated")
    gene_table = donor_table.groupby(["gene", "candidate", "run_role"], as_index=False).agg(
        n_donors=("subject", "nunique"),
        mean_overexpression_recovery_shift=("mean_overexpression_recovery_shift", "mean"),
        median_overexpression_recovery_shift=("mean_overexpression_recovery_shift", "median"),
        positive_donor_fraction=("mean_overexpression_recovery_shift", lambda v: float(np.mean(np.asarray(v) > 0))),
    )
    metrics_table = pd.DataFrame(metric_rows)
    gene_table = gene_table.merge(metrics_table, on=["gene", "candidate", "run_role"], how="left")
    gene_table.to_csv(output_dir / "gene_level_effects.tsv", sep="\t", index=False)

    gene_all_table = donor_all_table.groupby(["gene", "candidate", "run_role"], as_index=False).agg(
        n_donors=("subject", "nunique"),
        mean_overexpression_recovery_shift=("mean_overexpression_recovery_shift", "mean"),
        median_overexpression_recovery_shift=("mean_overexpression_recovery_shift", "median"),
        positive_donor_fraction=("mean_overexpression_recovery_shift", lambda v: float(np.mean(np.asarray(v) > 0))),
    )
    gene_all_table = gene_all_table.merge(metrics_table, on=["gene", "candidate", "run_role"], how="left")
    gene_all_table.to_csv(output_dir / "gene_level_effects_all_baseline.tsv", sep="\t", index=False)

    bootstrap_rows = []
    for gene, subset in donor_table.groupby("gene", sort=False):
        summary = bootstrap_summary(subset["mean_overexpression_recovery_shift"].to_numpy(float), rng)
        bootstrap_rows.append({"gene": gene, **summary})
    pd.DataFrame(bootstrap_rows).to_csv(output_dir / "bootstrap_summary.tsv", sep="\t", index=False)

    matched_rows = []
    for target in MVP_TARGETS:
        target_row = gene_table[(gene_table.gene == target) & (gene_table.run_role == "candidate")]
        controls = gene_table[(gene_table.candidate == target) & (gene_table.run_role == "matched_control")]
        if target_row.empty or len(controls) < 4:
            matched_rows.append({"candidate": target, "n_controls": len(controls), "candidate_effect": np.nan, "empirical_p_microtest": np.nan, "candidate_percentile_microtest": np.nan})
            continue
        effect = float(target_row.iloc[0]["mean_overexpression_recovery_shift"])
        null = controls["mean_overexpression_recovery_shift"].to_numpy(float)
        p = (1 + np.sum(null >= effect)) / (len(null) + 1)
        percentile = (np.sum(null < effect) + 0.5 * np.sum(null == effect)) / len(null)
        matched_rows.append({"candidate": target, "n_controls": len(null), "candidate_effect": effect, "control_mean": float(np.mean(null)), "control_sd": float(np.std(null, ddof=1)) if len(null) > 1 else np.nan, "empirical_p_microtest": float(p), "candidate_percentile_microtest": float(percentile)})
    pd.DataFrame(matched_rows).to_csv(output_dir / "matched_null_summary.tsv", sep="\t", index=False)

    # LODO rank stability across the small engineering set; not a formal biological inference.
    full_rank = donor_table.groupby("gene")["mean_overexpression_recovery_shift"].mean().sort_values(ascending=False)
    lodo_rows = []
    for donor in sorted(donor_table.subject.unique()):
        leave = donor_table[donor_table.subject != donor].groupby("gene")["mean_overexpression_recovery_shift"].mean()
        common = full_rank.index.intersection(leave.index)
        rho = spearmanr(full_rank.loc[common], leave.loc[common]).statistic if len(common) >= 3 else np.nan
        k = min(3, len(common))
        full_top = set(full_rank.loc[common].nlargest(k).index)
        leave_top = set(leave.loc[common].nlargest(k).index)
        jaccard = len(full_top & leave_top) / max(1, len(full_top | leave_top))
        lodo_rows.append({"left_out_subject": donor, "n_genes": len(common), "spearman_rank_stability": float(rho), "top3_jaccard": float(jaccard)})
    lodo_table = pd.DataFrame(lodo_rows)
    lodo_table.to_csv(output_dir / "lodo_rank_stability.tsv", sep="\t", index=False)

    baseline_by_donor, baseline_summary = simple_baselines(collapsed, meta, observed)
    baseline_by_donor.to_csv(output_dir / "simple_baseline_by_donor.tsv", sep="\t", index=False)
    baseline_summary.to_csv(output_dir / "simple_baselines.tsv", sep="\t", index=False)

    # The MVP gate is technical and stability-oriented. Matched-null significance is explicitly not a gate.
    observed_positive = float(np.mean(observed.geneformer_observed_recovery_shift > 0)) if len(observed) else 0.0
    median_lodo_rho = float(np.nanmedian(lodo_table.spearman_rank_stability))
    median_lodo_jaccard = float(np.nanmedian(lodo_table.top3_jaccard))
    candidate_donor_min = int(gene_table.loc[gene_table.run_role.eq("candidate"), "n_donors"].min())
    gate = {
        "decision_scope": "overexpression_engineering_MVP_only_not_formal_significance",
        "candidate_set_frozen": CANDIDATES,
        "mvp_targets": MVP_TARGETS,
        "n_balanced_cells": int(collapsed.shape[0]),
        "n_mapped_genes": int(collapsed.shape[1]),
        "n_lodo_donors": int(len(observed)),
        "observed_recovery_positive_fraction": observed_positive,
        "minimum_candidate_donor_coverage": candidate_donor_min,
        "median_lodo_rank_spearman": median_lodo_rho,
        "median_lodo_top3_jaccard": median_lodo_jaccard,
        "technical_pass": bool(len(observed) >= 12 and observed_positive >= 0.70 and candidate_donor_min >= 8 and np.isfinite(median_lodo_rho) and median_lodo_rho >= 0.60 and median_lodo_jaccard >= 0.50),
        "matched_null_significance_used_as_gate": False,
    }
    gate["next_step"] = "GO_TO_SEPARATELY_VERSIONED_FULL_10_PLUS_200" if gate["technical_pass"] else "STOP_AND_AUDIT_MVP"
    (output_dir / "mvp_gate.json").write_text(json.dumps(gate, indent=2, ensure_ascii=False), encoding="utf-8")

    qc = {
        "matrix_shape_input": list(x.shape),
        "matrix_shape_geneformer": list(collapsed.shape),
        "analysis_group_counts": {str(k): int(v) for k, v in group_counts.items()},
        "sequence_length_min": int(min(map(len, sequences))),
        "sequence_length_median": float(np.median(list(map(len, sequences)))),
        "sequence_length_max": int(max(map(len, sequences))),
        "model_class": model.__class__.__name__,
        "model_parameter_count": int(sum(p.numel() for p in model.parameters())),
        "device": str(device),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "all_embeddings_finite": bool(np.isfinite(embeddings).all()),
        "all_gene_effects_finite": bool(np.isfinite(donor_table.mean_overexpression_recovery_shift).all()),
        "all_baseline_gene_effects_finite": bool(np.isfinite(donor_all_table.mean_overexpression_recovery_shift).all()),
        "overexpression_sequence_contract": "V2 cls token retained first; target token moved or inserted immediately after cls; eos retained last; official overflow-matched original truncation applied",
        "model_max_input_length": model_max_len,
        "maximum_perturbed_sequence_length": maximum_perturbed_length,
        "n_overflow_matched_comparisons": int(overflow_matched_comparisons),
        "sequence_length_gate_pass": bool(maximum_perturbed_length <= model_max_len),
    }
    (output_dir / "technical_qc.json").write_text(json.dumps(qc, indent=2, ensure_ascii=False), encoding="utf-8")

    execution = {
        "seed": SEED,
        "bootstrap_reps": BOOTSTRAP_REPS,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "torch": torch.__version__,
        "transformers": __import__("transformers").__version__,
        "input_sha256": {str(p): sha256(p) for p in [matrix_path, genes_path, meta_path, manifest_path, MODEL_DIR / "model.safetensors", token_path, median_path, name_id_path]},
        "source_note": "One gene overexpressed at a time by moving or inserting its token immediately after V2 cls; no grouped perturbation; no endpoint/threshold/candidate tuning.",
        "perturbation": PERTURBATION,
    }
    (output_dir / "execution_manifest.json").write_text(json.dumps(execution, indent=2, ensure_ascii=False), encoding="utf-8")

    print("MVP_COMPLETE")
    print(json.dumps(gate, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
