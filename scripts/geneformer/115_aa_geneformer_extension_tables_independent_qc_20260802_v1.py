from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256(path.read_bytes())
    return digest.hexdigest()


def close(a: float, b: float, tolerance: float = 1e-12) -> bool:
    return bool(np.isclose(float(a), float(b), rtol=0, atol=tolerance))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables-dir", required=True)
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()
    global TABLES, OUT
    TABLES, OUT = Path(args.tables_dir), Path(args.output_json)
    a = pd.read_csv(TABLES / "Table_A_positive_control_bidirectional_perturbation.tsv", sep="\t")
    b = pd.read_csv(TABLES / "Table_B_state_specific_perturbation.tsv", sep="\t")
    c1 = pd.read_csv(TABLES / "Table_C1_simple_baselines.tsv", sep="\t")
    c2 = pd.read_csv(TABLES / "Table_C2_LODO_by_donor.tsv", sep="\t")
    c3 = pd.read_csv(TABLES / "Table_C3_LODO_summary.tsv", sep="\t")
    d = pd.read_csv(TABLES / "Table_D_program_ablation.tsv", sep="\t")
    e = pd.read_csv(TABLES / "Table_E_pretraining_leakage_audit.tsv", sep="\t")
    manifest = json.loads((TABLES / "manifest.json").read_text(encoding="utf-8"))

    checks: dict[str, bool] = {}
    checks["positive_controls_prespecified_six"] = list(a["gene"]) == [
        "MPL", "JAK2", "STAT5A", "STAT5B", "PIM1", "BCL2L1"
    ]
    checks["mpl_reason_matches_frozen_detection_audit"] = (
        a.loc[a["gene"].eq("MPL"), "measurement_status"].item() == "not_estimable"
        and "only_2_baseline_expressing_cells" in a.loc[a["gene"].eq("MPL"), "status_detail"].item()
    )
    checks["five_positive_controls_estimated"] = int(a["measurement_status"].eq("estimated").sum()) == 5
    checks["positive_matching_twenty_each"] = bool(
        (a.loc[a["measurement_status"].eq("estimated"), "n_matched_controls"] == 20).all()
    )
    checks["state_classes_exact"] = set(b["state_class"]) == {
        "HSPC-marker-class", "megakaryocyte-marker-class"
    }
    checks["megakaryocyte_scope_descriptive"] = bool(
        b.loc[b["state_class"].eq("megakaryocyte-marker-class"), "inference_scope"]
        .eq("descriptive_small_n_stress_test")
        .all()
    )
    checks["simple_baseline_rows_present"] = set(c1["baseline"]) == {
        "expression_centroid_shift",
        "ridge_shift",
        "geneformer_observed_recovery_shift",
        "spearman_expression_centroid_shift_vs_geneformer",
        "spearman_ridge_shift_vs_geneformer",
    }
    checks["lodo_summary_exact"] = (
        len(c2) == 15
        and close(c3["mean_spearman"].item(), c2["spearman_rank_stability"].mean())
        and close(c3["minimum_spearman"].item(), c2["spearman_rank_stability"].min())
        and close(c3["mean_top3_jaccard"].item(), c2["top3_jaccard"].mean())
    )
    checks["four_ablation_variants"] = set(d["variant"]) == {
        "all_programs",
        "leave_out_HSPC_identity",
        "leave_out_cell_cycle_recovery",
        "leave_out_hematopoietic_support",
    }
    checks["leakage_audit_has_residual_uncertainty"] = (
        "Leakage_conclusion" in set(e["audit_item"])
        and e.loc[e["audit_item"].eq("Sample_level_overlap"), "assessment"].item().startswith(
            "Exact sample-level overlap cannot be proven absent"
        )
    )
    checks["manifest_declares_no_scientific_input_modification"] = (
        manifest["scientific_inputs_modified"] is False
    )
    checks["manifest_hashes_match"] = all(
        sha256(TABLES / name) == digest for name, digest in manifest["output_sha256"].items()
    )

    result = {
        "all_checks_pass": all(checks.values()),
        "checks": checks,
        "verified_table_sha256": {
            path.name: sha256(path) for path in sorted(TABLES.glob("*.tsv"))
        },
    }
    if OUT.exists():
        raise RuntimeError(f"Refusing to overwrite existing QC artifact: {OUT}")
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if result["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
