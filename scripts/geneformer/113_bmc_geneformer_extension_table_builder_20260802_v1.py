from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(name: str) -> pd.DataFrame:
    return pd.read_csv(DATA / name, sep="\t")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    global DATA, OUT
    DATA, OUT = Path(args.input_dir), Path(args.output_dir)
    if OUT.exists() and any(OUT.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty output directory: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)

    pos_summary = read_tsv(
        "BMC_Geneformer_POSCTRL_MATCHED_20260802_v2/positive_control_gene_summary.tsv"
    )
    individual = read_tsv(
        "BMC_Geneformer_POSCTRL_MATCHED_20260802_v2/individual_matched_null.tsv"
    )
    individual = individual.rename(columns={"positive_control": "gene"})
    pos = pos_summary.merge(individual, on="gene", how="left", validate="one_to_one")
    pos["measurement_status"] = "estimated"
    pos["status_detail"] = "estimated_under_frozen_overall_gate"
    mpl = {column: pd.NA for column in pos.columns}
    mpl.update(
        {
            "gene": "MPL",
            "measurement_status": "not_estimable",
            "status_detail": (
                "only_2_baseline_expressing_cells; below_frozen_overall_"
                "expression_and_donor_coverage_gate"
            ),
        }
    )
    table_a = pd.concat([pos, pd.DataFrame([mpl])], ignore_index=True)
    order = {g: i for i, g in enumerate(["MPL", "JAK2", "STAT5A", "STAT5B", "PIM1", "BCL2L1"])}
    table_a["frozen_order"] = table_a["gene"].map(order)
    table_a = table_a.sort_values("frozen_order", kind="stable")
    table_a.to_csv(OUT / "Table_A_positive_control_bidirectional_perturbation.tsv", sep="\t", index=False)

    availability = read_tsv("BMC_Geneformer_STATE_SPECIFIC_20260802_v1_effect_availability.tsv")
    state_summary = read_tsv("BMC_Geneformer_STATE_SPECIFIC_20260802_v1_state_gene_summary.tsv")
    bootstrap = read_tsv("BMC_Geneformer_STATE_SPECIFIC_20260802_v1_state_gene_bootstrap.tsv")
    bootstrap = bootstrap.loc[bootstrap["metric"].eq("bidirectional_recovery_score")].copy()
    bootstrap = bootstrap.rename(
        columns={
            "mean": "bootstrap_mean_bidirectional_score",
            "median": "bootstrap_median_bidirectional_score",
            "ci_low": "bootstrap_ci_low",
            "ci_high": "bootstrap_ci_high",
            "positive_fraction": "bootstrap_positive_fraction",
        }
    )
    keep = [
        "state_class",
        "gene",
        "bootstrap_mean_bidirectional_score",
        "bootstrap_median_bidirectional_score",
        "bootstrap_ci_low",
        "bootstrap_ci_high",
        "bootstrap_positive_fraction",
    ]
    table_b = availability.merge(state_summary, on=["state_class", "gene", "n_donors"], how="left")
    table_b = table_b.merge(bootstrap[keep], on=["state_class", "gene"], how="left")
    table_b.to_csv(OUT / "Table_B_state_specific_perturbation.tsv", sep="\t", index=False)

    baseline = read_tsv("BMC_Geneformer_POSCTRL_DELETE_20260802_v2_simple_baselines.tsv")
    lodo = read_tsv("BMC_Geneformer_POSCTRL_DELETE_20260802_v2_lodo_rank_stability.tsv")
    lodo_summary = pd.DataFrame(
        [
            {
                "assessment": "LODO_gene_rank_stability",
                "n_left_out_donors": len(lodo),
                "mean_spearman": lodo["spearman_rank_stability"].mean(),
                "minimum_spearman": lodo["spearman_rank_stability"].min(),
                "mean_top3_jaccard": lodo["top3_jaccard"].mean(),
                "minimum_top3_jaccard": lodo["top3_jaccard"].min(),
            }
        ]
    )
    baseline.to_csv(OUT / "Table_C1_simple_baselines.tsv", sep="\t", index=False)
    lodo.to_csv(OUT / "Table_C2_LODO_by_donor.tsv", sep="\t", index=False)
    lodo_summary.to_csv(OUT / "Table_C3_LODO_summary.tsv", sep="\t", index=False)

    ablation = read_tsv("BMC_Geneformer_PROGRAM_ABLATION_20260802_v1_program_ablation_summary.tsv")
    ablation.to_csv(OUT / "Table_D_program_ablation.tsv", sep="\t", index=False)

    leakage = read_tsv("BMC_Geneformer_PRETRAINING_LEAKAGE_AUDIT_20260802_v1.tsv")
    leakage.to_csv(OUT / "Table_E_pretraining_leakage_audit.tsv", sep="\t", index=False)

    outputs = sorted(OUT.glob("*.tsv"))
    manifest = {
        "artifact_scope": "derived tables only; no figure or manuscript generated",
        "scientific_inputs_modified": False,
        "mpl_status_correction": (
            "MPL was present in the frozen positive-control list and detected in two baseline cells, "
            "but was not estimable under the frozen overall expression/donor-coverage gate."
        ),
        "output_sha256": {path.name: sha256(path) for path in outputs},
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
