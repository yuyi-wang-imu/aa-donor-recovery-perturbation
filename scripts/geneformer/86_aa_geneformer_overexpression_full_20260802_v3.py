from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ENGINE = Path(__file__).with_name("85_aa_geneformer_overexpression_mvp_20260802_v3.py")
spec = importlib.util.spec_from_file_location("aa_geneformer_overexpression_engine_v3", ENGINE)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Unable to load frozen Geneformer engine: {ENGINE}")
engine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(engine)


def find_preferred_source_asset(pattern: str) -> Path:
    preferred = engine.SOURCE_DIR / "geneformer" / pattern
    if preferred.exists():
        return preferred
    hits = [p for p in sorted(engine.SOURCE_DIR.rglob(pattern)) if "/build/lib/" not in p.as_posix()]
    if len(hits) != 1:
        raise RuntimeError(f"Expected one non-build source asset for {pattern}; observed {len(hits)}: {hits}")
    return hits[0]


def stratified_null(effects: pd.DataFrame, candidates: list[str], reps: int = 100000, seed: int = 20260802):
    rng = np.random.default_rng(seed)
    signed_candidates = []
    absolute_candidates = []
    control_vectors = []
    rows = []
    for candidate in candidates:
        cand = effects[(effects.gene == candidate) & (effects.run_role == "candidate")]
        ctrl = effects[(effects.candidate == candidate) & (effects.run_role == "matched_control")]
        if cand.empty or ctrl.empty:
            rows.append({"candidate": candidate, "status": "not_estimable", "n_controls": len(ctrl)})
            continue
        value = float(cand.iloc[0].mean_overexpression_recovery_shift)
        null = ctrl.mean_overexpression_recovery_shift.to_numpy(float)
        p_positive = (1 + np.sum(null >= value)) / (len(null) + 1)
        p_negative = (1 + np.sum(null <= value)) / (len(null) + 1)
        p_two_sided = min(1.0, 2 * min(p_positive, p_negative))
        p_absolute = (1 + np.sum(np.abs(null) >= abs(value))) / (len(null) + 1)
        rows.append({
            "candidate": candidate,
            "status": "estimated",
            "n_controls": len(null),
            "signed_effect": value,
            "absolute_effect": abs(value),
            "control_mean": float(np.mean(null)),
            "control_sd": float(np.std(null, ddof=1)) if len(null) > 1 else np.nan,
            "p_toward_healthy": float(p_positive),
            "p_away_from_healthy": float(p_negative),
            "p_two_sided": float(p_two_sided),
            "p_absolute_magnitude": float(p_absolute),
        })
        signed_candidates.append(value)
        absolute_candidates.append(abs(value))
        control_vectors.append(null)
    if not control_vectors:
        return pd.DataFrame(rows), {"status": "not_estimable"}
    observed_signed = float(np.mean(signed_candidates))
    observed_absolute = float(np.mean(absolute_candidates))
    null_signed = np.empty(reps)
    null_absolute = np.empty(reps)
    for i in range(reps):
        sampled = np.asarray([rng.choice(vector) for vector in control_vectors], dtype=float)
        null_signed[i] = np.mean(sampled)
        null_absolute[i] = np.mean(np.abs(sampled))
    pooled = {
        "status": "estimated",
        "n_candidate_strata": len(control_vectors),
        "null_reps": reps,
        "observed_mean_signed_effect": observed_signed,
        "observed_mean_absolute_effect": observed_absolute,
        "empirical_p_toward_healthy": float((1 + np.sum(null_signed >= observed_signed)) / (reps + 1)),
        "empirical_p_away_from_healthy": float((1 + np.sum(null_signed <= observed_signed)) / (reps + 1)),
        "empirical_p_two_sided_signed": float(min(1.0, 2 * min((1 + np.sum(null_signed >= observed_signed)) / (reps + 1), (1 + np.sum(null_signed <= observed_signed)) / (reps + 1)))),
        "empirical_p_absolute_magnitude": float((1 + np.sum(null_absolute >= observed_absolute)) / (reps + 1)),
        "null_signed_q025": float(np.quantile(null_signed, 0.025)),
        "null_signed_q975": float(np.quantile(null_signed, 0.975)),
        "null_absolute_q025": float(np.quantile(null_absolute, 0.025)),
        "null_absolute_q975": float(np.quantile(null_absolute, 0.975)),
    }
    return pd.DataFrame(rows), pooled


def matching_diagnostics(effects: pd.DataFrame, candidates: list[str]):
    covariates = ["baseline_detection_fraction", "median_normalized_token_rank"]
    rows = []
    for candidate in candidates:
        cand = effects[(effects.gene == candidate) & (effects.run_role == "candidate")]
        ctrl = effects[(effects.candidate == candidate) & (effects.run_role == "matched_control")]
        if cand.empty or ctrl.empty:
            rows.append({"candidate": candidate, "status": "not_estimable"})
            continue
        record = {"candidate": candidate, "status": "estimated", "n_controls": len(ctrl)}
        imbalances = []
        for covariate in covariates:
            target = float(cand.iloc[0][covariate])
            values = ctrl[covariate].to_numpy(float)
            sd = float(np.std(values, ddof=1)) if len(values) > 1 else np.nan
            residual = abs(target - float(np.mean(values))) / sd if np.isfinite(sd) and sd > 0 else np.nan
            record[f"{covariate}_target"] = target
            record[f"{covariate}_control_mean"] = float(np.mean(values))
            record[f"{covariate}_standardized_residual"] = residual
            if np.isfinite(residual): imbalances.append(residual)
        record["max_abs_standardized_residual"] = max(imbalances) if imbalances else np.nan
        record["formal_matched_claim_allowed"] = bool(imbalances and max(imbalances) <= 1.0)
        rows.append(record)
    return pd.DataFrame(rows)


def sctenifold_consistency(effects: pd.DataFrame, candidates: list[str]):
    endpoint_root = engine.CONTROL_ROOT / "runs/candidate"
    records = []
    for candidate in candidates:
        cand = effects[(effects.gene == candidate) & (effects.run_role == "candidate")]
        hits = sorted(endpoint_root.glob(f"*_candidate_{candidate}_endpoint.csv"))
        if cand.empty or len(hits) != 1:
            continue
        endpoint = pd.read_csv(hits[0]).iloc[0]
        records.append({
            "candidate": candidate,
            "geneformer_signed_overexpression_shift": float(cand.iloc[0].mean_overexpression_recovery_shift),
            "geneformer_absolute_overexpression_shift": abs(float(cand.iloc[0].mean_overexpression_recovery_shift)),
            "sctenifold_significant_response_genes": float(endpoint.n_sig_excluding_gKO_padj_0_05),
            "sctenifold_max_abs_Z": float(endpoint.max_abs_Z_excluding_gKO_exploratory),
        })
    table = pd.DataFrame(records)
    summary = {}
    if len(table) >= 3:
        summary = {
            "n_candidates": len(table),
            "spearman_signed_vs_response_count": float(spearmanr(table.geneformer_signed_overexpression_shift, table.sctenifold_significant_response_genes).statistic),
            "spearman_absolute_vs_response_count": float(spearmanr(table.geneformer_absolute_overexpression_shift, table.sctenifold_significant_response_genes).statistic),
            "spearman_absolute_vs_max_abs_Z": float(spearmanr(table.geneformer_absolute_overexpression_shift, table.sctenifold_max_abs_Z).statistic),
        }
    return table, summary


def main() -> int:
    # The audited engine is reused without editing.  Full scope is injected before execution.
    engine.find_one = find_preferred_source_asset
    engine.MVP_TARGETS = list(engine.CANDIDATES)
    engine.MAX_CONTROLS_PER_TARGET = 20
    original_argv = list(sys.argv)
    sys.argv = [original_argv[0]] + ["mvp" if arg == "full" else arg for arg in original_argv[1:]]
    rc = engine.main()
    sys.argv = original_argv
    if rc != 0:
        return rc

    out_arg = original_argv[original_argv.index("--output-dir") + 1]
    output_dir = Path(out_arg)
    effects = pd.read_csv(output_dir / "gene_level_effects.tsv", sep="\t")
    candidate_null, pooled = stratified_null(effects, engine.CANDIDATES)
    candidate_null.to_csv(output_dir / "candidate_stratified_matched_null.tsv", sep="\t", index=False)
    (output_dir / "pooled_stratified_matched_null.json").write_text(json.dumps(pooled, indent=2), encoding="utf-8")

    matching = matching_diagnostics(effects, engine.CANDIDATES)
    matching.to_csv(output_dir / "geneformer_matching_diagnostics.tsv", sep="\t", index=False)

    consistency, consistency_summary = sctenifold_consistency(effects, engine.CANDIDATES)
    consistency.to_csv(output_dir / "geneformer_sctenifold_consistency.tsv", sep="\t", index=False)
    (output_dir / "geneformer_sctenifold_consistency_summary.json").write_text(json.dumps(consistency_summary, indent=2), encoding="utf-8")

    estimated = candidate_null[candidate_null.status.eq("estimated")]
    controls_total = int(effects.loc[effects.run_role.eq("matched_control"), "gene"].nunique())
    full_gate = {
        "decision_scope": "overexpression_full_10_candidate_plus_frozen_matched_controls",
        "candidate_set_frozen": engine.CANDIDATES,
        "n_candidates_estimated": int(len(estimated)),
        "n_unique_controls_estimated": controls_total,
        "n_candidate_strata_with_20_controls": int((candidate_null.n_controls.fillna(0) == 20).sum()),
        "pooled_matched_null": pooled,
        "matching_diagnostic_only": True,
        "formal_matched_claim_allowed_for_n_candidates": int(matching.formal_matched_claim_allowed.fillna(False).sum()),
        "technical_pass": bool(len(estimated) >= 9 and controls_total >= 180),
        "significance_used_as_technical_gate": False,
    }
    full_gate["next_step"] = "FREEZE_OVEREXPRESSION_RESULTS_AND_RUN_BIDIRECTIONAL_AUDIT" if full_gate["technical_pass"] else "STOP_AND_AUDIT_OVEREXPRESSION_COVERAGE"
    (output_dir / "full_analysis_gate.json").write_text(json.dumps(full_gate, indent=2), encoding="utf-8")
    print("OVEREXPRESSION_FULL_ANALYSIS_COMPLETE")
    print(json.dumps(full_gate))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
