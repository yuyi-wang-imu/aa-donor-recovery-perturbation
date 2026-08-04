from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 20260802
K_CONTROLS = 20
POOL_REPS = 100000
BALANCE_LIMIT = 1.0
CANDIDATES = ["CDK6", "CA2", "PARP1", "KIT", "SYK", "GSK3B", "HIF1A", "TOP2A", "TERT", "CD38"]
COVARIATES = ["baseline_detection_fraction", "median_normalized_token_rank"]


def refuse_existing(paths: list[Path]) -> None:
    hits = [path for path in paths if path.exists()]
    if hits:
        raise FileExistsError(f"Refusing to overwrite existing outputs: {hits}")


def bh_adjust(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(float)
    result = np.full(numeric.shape, np.nan, dtype=float)
    finite = np.flatnonzero(np.isfinite(numeric))
    if not len(finite):
        return pd.Series(result, index=values.index)
    order = finite[np.argsort(numeric[finite], kind="mergesort")]
    m = len(order)
    adjusted = np.minimum.accumulate((numeric[order] * m / np.arange(1, m + 1))[::-1])[::-1]
    result[order] = np.minimum(adjusted, 1.0)
    return pd.Series(result, index=values.index)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--effects", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    effects_path = Path(args.effects)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = "BMC_Geneformer_REMATCH_20260802_v1"
    outputs = {
        "selection": output_dir / f"{stem}_selection.tsv",
        "diagnostics": output_dir / f"{stem}_diagnostics.tsv",
        "candidate_null": output_dir / f"{stem}_candidate_matched_null.tsv",
        "pooled": output_dir / f"{stem}_pooled_matched_null.json",
        "gate": output_dir / f"{stem}_gate.json",
    }
    refuse_existing(list(outputs.values()))

    effects = pd.read_csv(effects_path, sep="\t")
    required = {"gene", "candidate", "run_role", "mean_deletion_recovery_shift", *COVARIATES}
    missing = required.difference(effects.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    effects[COVARIATES] = effects[COVARIATES].apply(pd.to_numeric, errors="raise")
    candidates = effects.loc[effects.run_role.eq("candidate")].drop_duplicates("gene").set_index("gene")
    pool = effects.loc[effects.run_role.eq("matched_control")].drop_duplicates("gene").copy()
    pool = pool.loc[~pool.gene.isin(CANDIDATES)].reset_index(drop=True)
    if not set(CANDIDATES).issubset(candidates.index):
        raise ValueError(f"Missing frozen candidates: {sorted(set(CANDIDATES).difference(candidates.index))}")
    if len(pool) < K_CONTROLS:
        raise ValueError(f"Insufficient unique controls: {len(pool)}")

    scale = pool[COVARIATES].std(ddof=1).to_numpy(float)
    if np.any(~np.isfinite(scale)) or np.any(scale <= 0):
        raise ValueError(f"Invalid covariate scale: {scale}")

    selection_rows: list[dict] = []
    diagnostic_rows: list[dict] = []
    null_rows: list[dict] = []
    eligible_vectors: list[np.ndarray] = []
    eligible_candidate_effects: list[float] = []

    for candidate in CANDIDATES:
        target = candidates.loc[candidate]
        target_cov = target[COVARIATES].to_numpy(float)
        delta = (pool[COVARIATES].to_numpy(float) - target_cov) / scale
        distance = np.sqrt(np.sum(delta * delta, axis=1))
        selected = (
            pool.assign(match_distance=distance)
            .sort_values(["match_distance", "gene"], kind="mergesort")
            .head(K_CONTROLS)
            .copy()
        )
        for rank, row in enumerate(selected.itertuples(index=False), start=1):
            selection_rows.append({
                "candidate": candidate,
                "control_gene": row.gene,
                "match_rank": rank,
                "match_distance": float(row.match_distance),
                "baseline_detection_fraction": float(row.baseline_detection_fraction),
                "median_normalized_token_rank": float(row.median_normalized_token_rank),
                "control_effect": float(row.mean_deletion_recovery_shift),
            })

        diagnostic = {"candidate": candidate, "n_controls": len(selected)}
        residuals = []
        for covariate in COVARIATES:
            target_value = float(target[covariate])
            control_mean = float(selected[covariate].mean())
            control_sd = float(selected[covariate].std(ddof=1))
            residual = abs(target_value - control_mean) / control_sd if control_sd > 0 else np.nan
            diagnostic[f"{covariate}_target"] = target_value
            diagnostic[f"{covariate}_control_mean"] = control_mean
            diagnostic[f"{covariate}_standardized_residual"] = residual
            if np.isfinite(residual):
                residuals.append(residual)
        diagnostic["max_abs_standardized_residual"] = max(residuals) if residuals else np.nan
        diagnostic["formal_matched_claim_allowed"] = bool(
            residuals and diagnostic["max_abs_standardized_residual"] <= BALANCE_LIMIT
        )
        diagnostic_rows.append(diagnostic)

        candidate_effect = float(target.mean_deletion_recovery_shift)
        null = selected.mean_deletion_recovery_shift.to_numpy(float)
        p_toward = float((1 + np.sum(null >= candidate_effect)) / (len(null) + 1))
        p_away = float((1 + np.sum(null <= candidate_effect)) / (len(null) + 1))
        record = {
            "candidate": candidate,
            "n_controls": len(null),
            "candidate_effect": candidate_effect,
            "control_mean": float(np.mean(null)),
            "control_sd": float(np.std(null, ddof=1)),
            "p_toward_healthy": p_toward,
            "p_away_from_healthy": p_away,
            "p_two_sided": min(1.0, 2 * min(p_toward, p_away)),
            "p_absolute_magnitude": float((1 + np.sum(np.abs(null) >= abs(candidate_effect))) / (len(null) + 1)),
            "formal_matched_claim_allowed": diagnostic["formal_matched_claim_allowed"],
        }
        null_rows.append(record)
        if diagnostic["formal_matched_claim_allowed"]:
            eligible_vectors.append(null)
            eligible_candidate_effects.append(candidate_effect)

    selection = pd.DataFrame(selection_rows)
    diagnostics = pd.DataFrame(diagnostic_rows)
    candidate_null = pd.DataFrame(null_rows)
    for column in ["p_toward_healthy", "p_away_from_healthy", "p_two_sided", "p_absolute_magnitude"]:
        candidate_null[f"{column}_bh"] = bh_adjust(candidate_null[column])

    rng = np.random.default_rng(SEED)
    pooled = {
        "status": "not_estimable",
        "matching_uses_outcomes": False,
        "matching_covariates": COVARIATES,
        "k_controls_per_candidate": K_CONTROLS,
        "balance_limit": BALANCE_LIMIT,
        "null_reps": POOL_REPS,
        "n_eligible_candidate_strata": len(eligible_vectors),
    }
    if eligible_vectors:
        observed_signed = float(np.mean(eligible_candidate_effects))
        observed_absolute = float(np.mean(np.abs(eligible_candidate_effects)))
        null_signed = np.empty(POOL_REPS, dtype=float)
        null_absolute = np.empty(POOL_REPS, dtype=float)
        for i in range(POOL_REPS):
            sampled = np.asarray([rng.choice(vector) for vector in eligible_vectors], dtype=float)
            null_signed[i] = np.mean(sampled)
            null_absolute[i] = np.mean(np.abs(sampled))
        p_toward = float((1 + np.sum(null_signed >= observed_signed)) / (POOL_REPS + 1))
        p_away = float((1 + np.sum(null_signed <= observed_signed)) / (POOL_REPS + 1))
        pooled.update({
            "status": "estimated",
            "observed_mean_signed_effect": observed_signed,
            "observed_mean_absolute_effect": observed_absolute,
            "empirical_p_toward_healthy": p_toward,
            "empirical_p_away_from_healthy": p_away,
            "empirical_p_two_sided_signed": min(1.0, 2 * min(p_toward, p_away)),
            "empirical_p_absolute_magnitude": float((1 + np.sum(null_absolute >= observed_absolute)) / (POOL_REPS + 1)),
            "null_signed_q025": float(np.quantile(null_signed, 0.025)),
            "null_signed_q975": float(np.quantile(null_signed, 0.975)),
            "null_absolute_q025": float(np.quantile(null_absolute, 0.025)),
            "null_absolute_q975": float(np.quantile(null_absolute, 0.975)),
        })

    gate = {
        "scope": "outcome_blind_covariate_rematch_within_already_computed_control_union",
        "candidate_set_frozen": CANDIDATES,
        "n_candidate_strata": len(CANDIDATES),
        "n_unique_available_controls": int(pool.gene.nunique()),
        "n_unique_selected_controls": int(selection.control_gene.nunique()),
        "n_formally_balanced_strata": int(diagnostics.formal_matched_claim_allowed.sum()),
        "unbalanced_strata": diagnostics.loc[~diagnostics.formal_matched_claim_allowed, "candidate"].tolist(),
        "technical_pass": bool(diagnostics.formal_matched_claim_allowed.sum() >= 8),
        "significance_used_for_matching_or_gate": False,
        "pooled_analysis_restricted_to_balanced_strata": True,
    }

    selection.to_csv(outputs["selection"], sep="\t", index=False)
    diagnostics.to_csv(outputs["diagnostics"], sep="\t", index=False)
    candidate_null.to_csv(outputs["candidate_null"], sep="\t", index=False)
    outputs["pooled"].write_text(json.dumps(pooled, indent=2), encoding="utf-8")
    outputs["gate"].write_text(json.dumps(gate, indent=2), encoding="utf-8")
    print(json.dumps(gate))
    print(json.dumps(pooled))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
