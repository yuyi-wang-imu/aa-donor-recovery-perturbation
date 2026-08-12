from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 20260802
BOOTSTRAP_REPS = 2000
POOLED_NULL_REPS = 100000
SPEARMAN_PERMUTATION_REPS = 100000
CANDIDATES = ["CDK6", "CA2", "PARP1", "KIT", "SYK", "GSK3B", "HIF1A", "TOP2A", "TERT", "CD38"]


def refuse_nonempty(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise RuntimeError(f"Refusing to write into non-empty output directory: {path}")
    path.mkdir(parents=True, exist_ok=True)


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(1.0, adjusted)
    return output


def collapse_gene_effects(table: pd.DataFrame, column: str) -> pd.Series:
    spread = table.groupby("gene")[column].agg(lambda values: float(np.max(values) - np.min(values)))
    if (spread > 1e-10).any():
        raise RuntimeError(f"Gene effects vary across duplicated matching strata: {spread[spread > 1e-10].to_dict()}")
    return table.groupby("gene")[column].first()


def spearman_with_fixed_permutation(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    x_rank = pd.Series(x).rank(method="average").to_numpy(float)
    y_rank = pd.Series(y).rank(method="average").to_numpy(float)
    x_centered = x_rank - x_rank.mean()
    y_centered = y_rank - y_rank.mean()
    denominator = float(np.sqrt(np.sum(x_centered**2) * np.sum(y_centered**2)))
    if denominator == 0:
        return float("nan"), float("nan")
    observed = float(np.dot(x_centered, y_centered) / denominator)
    rng = np.random.default_rng(SEED + 17)
    exceed = 0
    completed = 0
    batch_size = 5000
    while completed < SPEARMAN_PERMUTATION_REPS:
        size = min(batch_size, SPEARMAN_PERMUTATION_REPS - completed)
        order = np.argsort(rng.random((size, len(y_centered))), axis=1)
        permuted = y_centered[order]
        correlations = (permuted @ x_centered) / denominator
        exceed += int(np.sum(np.abs(correlations) >= abs(observed) - 1e-15))
        completed += size
    return observed, float((1 + exceed) / (SPEARMAN_PERMUTATION_REPS + 1))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--deletion-gene-effects", required=True)
    parser.add_argument("--deletion-donor-effects", required=True)
    parser.add_argument("--overexpression-gene-effects", required=True)
    parser.add_argument("--overexpression-donor-effects", required=True)
    parser.add_argument("--overexpression-all-baseline-effects", required=True)
    parser.add_argument("--rematch-selection", required=True)
    parser.add_argument("--rematch-diagnostics", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    refuse_nonempty(output_dir)
    deletion = pd.read_csv(args.deletion_gene_effects, sep="\t")
    deletion_donor = pd.read_csv(args.deletion_donor_effects, sep="\t")
    overexpression = pd.read_csv(args.overexpression_gene_effects, sep="\t")
    overexpression_donor = pd.read_csv(args.overexpression_donor_effects, sep="\t")
    overexpression_all = pd.read_csv(args.overexpression_all_baseline_effects, sep="\t")
    rematch = pd.read_csv(args.rematch_selection, sep="\t")
    diagnostics = pd.read_csv(args.rematch_diagnostics, sep="\t")

    del_gene = collapse_gene_effects(deletion, "mean_deletion_recovery_shift")
    oe_gene = collapse_gene_effects(overexpression, "mean_overexpression_recovery_shift")
    oe_all_gene = collapse_gene_effects(overexpression_all, "mean_overexpression_recovery_shift")
    missing = sorted((set(CANDIDATES) - set(del_gene.index)) | (set(CANDIDATES) - set(oe_gene.index)))
    if missing:
        raise RuntimeError(f"Missing frozen candidates: {missing}")

    candidate_rows = []
    null_vectors: dict[str, np.ndarray] = {}
    for candidate in CANDIDATES:
        deletion_effect = float(del_gene[candidate])
        overexpression_effect = float(oe_gene[candidate])
        score = min(overexpression_effect, -deletion_effect)
        control_genes = rematch.loc[rematch.candidate.eq(candidate), "control_gene"].astype(str).tolist()
        unavailable = sorted((set(control_genes) - set(del_gene.index)) | (set(control_genes) - set(oe_gene.index)))
        if unavailable:
            raise RuntimeError(f"Missing frozen controls for {candidate}: {unavailable}")
        control_scores = np.asarray([min(float(oe_gene[gene]), -float(del_gene[gene])) for gene in control_genes])
        null_vectors[candidate] = control_scores
        p_enrichment = float((1 + np.sum(control_scores >= score)) / (len(control_scores) + 1))
        balance = diagnostics.loc[diagnostics.candidate.eq(candidate)].iloc[0]
        formal_allowed = str(balance.formal_matched_claim_allowed).strip().lower() == "true"
        candidate_rows.append({
            "candidate": candidate,
            "deletion_recovery_shift": deletion_effect,
            "overexpression_recovery_shift_detected_cells": overexpression_effect,
            "overexpression_recovery_shift_all_baseline_cells": float(oe_all_gene[candidate]),
            "bidirectional_score_min_of_arms": score,
            "directionally_coherent": bool(deletion_effect < 0 and overexpression_effect > 0),
            "contrast_overexpression_minus_deletion": overexpression_effect - deletion_effect,
            "n_frozen_matched_controls": len(control_scores),
            "matched_control_score_mean": float(np.mean(control_scores)),
            "matched_control_score_sd": float(np.std(control_scores, ddof=1)),
            "empirical_p_bidirectional_enrichment": p_enrichment,
            "matched_control_percentile": float(np.mean(control_scores <= score)),
            "max_abs_standardized_residual": float(balance.max_abs_standardized_residual),
            "formal_matched_claim_allowed": formal_allowed,
        })
    candidate_table = pd.DataFrame(candidate_rows)
    candidate_table["bh_q_bidirectional_enrichment"] = bh_adjust(candidate_table.empirical_p_bidirectional_enrichment.to_numpy(float))
    candidate_table.to_csv(output_dir / "candidate_bidirectional_summary.tsv", sep="\t", index=False)

    del_donor = deletion_donor[["gene", "subject", "mean_deletion_recovery_shift"]].drop_duplicates()
    oe_donor = overexpression_donor[["gene", "subject", "mean_overexpression_recovery_shift"]].drop_duplicates()
    donor = del_donor.merge(oe_donor, on=["gene", "subject"], how="inner")
    donor = donor[donor.gene.isin(CANDIDATES)].copy()
    donor["bidirectional_score_min_of_arms"] = np.minimum(donor.mean_overexpression_recovery_shift, -donor.mean_deletion_recovery_shift)
    donor["directionally_coherent"] = (donor.mean_deletion_recovery_shift < 0) & (donor.mean_overexpression_recovery_shift > 0)
    donor["contrast_overexpression_minus_deletion"] = donor.mean_overexpression_recovery_shift - donor.mean_deletion_recovery_shift
    donor.to_csv(output_dir / "donor_paired_bidirectional_effects.tsv", sep="\t", index=False)

    rng = np.random.default_rng(SEED)
    bootstrap_rows = []
    for candidate in CANDIDATES:
        values = donor.loc[donor.gene.eq(candidate), "bidirectional_score_min_of_arms"].to_numpy(float)
        boot = np.asarray([np.mean(rng.choice(values, size=len(values), replace=True)) for _ in range(BOOTSTRAP_REPS)])
        bootstrap_rows.append({
            "candidate": candidate,
            "n_paired_donors": len(values),
            "mean_donor_bidirectional_score": float(np.mean(values)),
            "bootstrap_ci025": float(np.quantile(boot, 0.025)),
            "bootstrap_ci975": float(np.quantile(boot, 0.975)),
            "donor_directional_coherence_fraction": float(donor.loc[donor.gene.eq(candidate), "directionally_coherent"].mean()),
        })
    pd.DataFrame(bootstrap_rows).to_csv(output_dir / "candidate_donor_bootstrap.tsv", sep="\t", index=False)

    balanced = candidate_table.loc[candidate_table.formal_matched_claim_allowed, "candidate"].tolist()
    observed = float(candidate_table.loc[candidate_table.candidate.isin(balanced), "bidirectional_score_min_of_arms"].mean())
    pooled_null = np.asarray([
        np.mean([rng.choice(null_vectors[candidate]) for candidate in balanced])
        for _ in range(POOLED_NULL_REPS)
    ])
    pooled = {
        "candidate_set_frozen": CANDIDATES,
        "balanced_candidate_strata": balanced,
        "unbalanced_strata_abstained": sorted(set(CANDIDATES) - set(balanced)),
        "bidirectional_score_definition": "min(overexpression recovery shift, negative deletion recovery shift)",
        "observed_mean_bidirectional_score": observed,
        "null_reps": POOLED_NULL_REPS,
        "empirical_p_enrichment": float((1 + np.sum(pooled_null >= observed)) / (POOLED_NULL_REPS + 1)),
        "null_q025": float(np.quantile(pooled_null, 0.025)),
        "null_q975": float(np.quantile(pooled_null, 0.975)),
    }
    (output_dir / "pooled_bidirectional_matched_null.json").write_text(json.dumps(pooled, indent=2), encoding="utf-8")

    correlation, correlation_p = spearman_with_fixed_permutation(
        candidate_table.overexpression_recovery_shift_detected_cells.to_numpy(float),
        -candidate_table.deletion_recovery_shift.to_numpy(float),
    )
    gate = {
        "scientific_question": "Do frozen candidates show donor-stable bidirectional consistency along the hematopoietic recovery axis beyond matched background?",
        "n_candidates": len(candidate_table),
        "n_directionally_coherent_candidates": int(candidate_table.directionally_coherent.sum()),
        "n_balanced_candidate_strata": len(balanced),
        "spearman_overexpression_vs_negative_deletion": correlation,
        "spearman_p_value_fixed_100000_permutations_descriptive": correlation_p,
        "pooled_matched_null": pooled,
        "technical_pass": bool(len(candidate_table) == 10 and len(balanced) == 8 and np.isfinite(candidate_table.bidirectional_score_min_of_arms).all()),
        "significance_used_as_technical_gate": False,
        "candidate_or_threshold_tuning_after_results": False,
        "scipy_dependency_required": False,
    }
    (output_dir / "bidirectional_audit_gate.json").write_text(json.dumps(gate, indent=2), encoding="utf-8")
    print(json.dumps(gate))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
