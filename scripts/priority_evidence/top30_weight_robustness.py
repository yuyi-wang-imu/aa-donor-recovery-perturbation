#!/usr/bin/env python3
"""Evaluate the stability of the frozen 30-candidate list across recorded weights.

This script extends the previously archived non-docking sensitivity analysis from
Top20 to Top30. It does not redefine the formal score or change the frozen order.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


COMPONENTS = [
    "score_WGCNA_module_strength",
    "score_module_pathology_importance",
    "score_CD34_expression",
    "score_CD34_marker_context",
    "score_AA_direction_recovery",
    "score_GSE165870_bulk_support",
    "score_network_pharmacology_hit",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", required=True, type=Path)
    parser.add_argument("--scenarios", required=True, type=Path)
    parser.add_argument("--formal-top30", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def is_true(value: object) -> bool:
    return str(value).strip().upper() == "TRUE"


def main() -> None:
    args = parse_args()
    scores = pd.read_csv(args.scores)
    scenarios = pd.read_csv(args.scenarios)
    formal = pd.read_csv(args.formal_top30)

    missing = sorted(set(COMPONENTS) - set(scores.columns))
    if missing:
        raise ValueError(f"Missing score columns: {missing}")

    eligible = scores[
        scores["module_robust_after_tech_adjustment"].map(is_true)
        & scores["technical_confounding_flag"].fillna("").ne("TECH_CAUTION")
    ].copy()
    if eligible.empty:
        raise ValueError("No eligible candidates remain after the archived filters.")

    formal_order = formal["GeneSymbol"].astype(str).tolist()
    if len(formal_order) != 30 or len(set(formal_order)) != 30:
        raise ValueError("The formal candidate file must contain 30 unique genes.")
    formal_rank = {gene: rank for rank, gene in enumerate(formal_order, start=1)}

    scenario_rows: list[dict[str, object]] = []
    selected_records: list[dict[str, object]] = []

    for row in scenarios.itertuples(index=False):
        weights = json.loads(row.weights_json)
        # Leave-one-component-out scenarios omit the removed component from JSON.
        weights = {column: float(weights.get(column, 0.0)) for column in COMPONENTS}

        candidate = eligible.copy()
        candidate["scenario_score"] = sum(
            pd.to_numeric(candidate[column], errors="coerce").fillna(0.0)
            * float(weights[column])
            for column in COMPONENTS
        )
        candidate = candidate.sort_values(
            ["scenario_score", "experimental_evidence_score", "GeneSymbol"],
            ascending=[False, False, True],
            kind="mergesort",
        ).reset_index(drop=True)
        candidate["scenario_rank"] = np.arange(1, len(candidate) + 1)
        top30 = candidate.head(30).copy()
        selected = top30["GeneSymbol"].astype(str).tolist()
        overlap_genes = [gene for gene in selected if gene in formal_rank]
        overlap = len(overlap_genes)
        union = len(set(selected).union(formal_order))
        jaccard = overlap / union if union else np.nan

        rank_pairs = pd.DataFrame(
            {
                "formal_rank": [formal_rank[gene] for gene in overlap_genes],
                "scenario_rank": [
                    int(
                        top30.loc[
                            top30["GeneSymbol"].astype(str).eq(gene), "scenario_rank"
                        ].iloc[0]
                    )
                    for gene in overlap_genes
                ],
            }
        )
        spearman = (
            rank_pairs["formal_rank"]
            .rank(method="average")
            .corr(rank_pairs["scenario_rank"].rank(method="average"))
            if len(rank_pairs) >= 3
            else np.nan
        )

        scenario_rows.append(
            {
                "scenario_id": row.scenario_id,
                "scenario_type": row.scenario_type,
                "note": row.note,
                "top30_overlap_count": overlap,
                "top30_overlap_fraction": overlap / 30.0,
                "top30_jaccard": jaccard,
                "spearman_rank_within_overlap": spearman,
                "top30_genes": ", ".join(selected),
                "weights_json": row.weights_json,
            }
        )
        for item in top30.itertuples(index=False):
            selected_records.append(
                {
                    "scenario_id": row.scenario_id,
                    "scenario_type": row.scenario_type,
                    "GeneSymbol": item.GeneSymbol,
                    "scenario_rank": int(item.scenario_rank),
                    "scenario_score": float(item.scenario_score),
                    "formal_rank": formal_rank.get(str(item.GeneSymbol), np.nan),
                    "in_formal_top30": str(item.GeneSymbol) in formal_rank,
                }
            )

    scenario_out = pd.DataFrame(scenario_rows)
    selected_out = pd.DataFrame(selected_records)

    all_genes = scores[["GeneSymbol", "rank_no_docking_formal"]].copy()
    freq = (
        selected_out.groupby("GeneSymbol", as_index=False)
        .agg(
            selected_scenarios=("scenario_id", "nunique"),
            median_selected_rank=("scenario_rank", "median"),
            minimum_selected_rank=("scenario_rank", "min"),
            maximum_selected_rank=("scenario_rank", "max"),
        )
    )
    gene_out = all_genes.merge(freq, on="GeneSymbol", how="left")
    gene_out["selected_scenarios"] = gene_out["selected_scenarios"].fillna(0).astype(int)
    gene_out["selection_frequency"] = (
        gene_out["selected_scenarios"] / scenario_out["scenario_id"].nunique()
    )
    gene_out["in_formal_top30"] = gene_out["GeneSymbol"].isin(formal_order)
    gene_out = gene_out.sort_values(
        ["in_formal_top30", "rank_no_docking_formal"],
        ascending=[False, True],
        kind="mergesort",
    )

    def summarize(frame: pd.DataFrame, label: str) -> dict[str, object]:
        return {
            "scenario_type": label,
            "n_scenarios": len(frame),
            "minimum_overlap_count": int(frame["top30_overlap_count"].min()),
            "median_overlap_count": float(frame["top30_overlap_count"].median()),
            "mean_overlap_count": float(frame["top30_overlap_count"].mean()),
            "maximum_overlap_count": int(frame["top30_overlap_count"].max()),
            "minimum_jaccard": float(frame["top30_jaccard"].min()),
            "median_jaccard": float(frame["top30_jaccard"].median()),
            "minimum_spearman_within_overlap": float(
                frame["spearman_rank_within_overlap"].min()
            ),
            "median_spearman_within_overlap": float(
                frame["spearman_rank_within_overlap"].median()
            ),
        }

    summary_rows = [summarize(scenario_out, "all_scenarios")]
    for scenario_type, frame in scenario_out.groupby("scenario_type", sort=True):
        summary_rows.append(summarize(frame, str(scenario_type)))
    summary_out = pd.DataFrame(summary_rows)

    derivation = pd.DataFrame(
        [
            (
                "WGCNA module strength",
                "score_WGCNA_module_strength",
                0.30,
                0.23478260869565215,
            ),
            (
                "Module pathology importance",
                "score_module_pathology_importance",
                0.18,
                0.14086956521739130,
            ),
            (
                "CD34 expression",
                "score_CD34_expression",
                0.18,
                0.14086956521739130,
            ),
            (
                "CD34 marker context",
                "score_CD34_marker_context",
                0.14,
                0.10956521739130434,
            ),
            (
                "AA-direction recovery",
                "score_AA_direction_recovery",
                0.14,
                0.10956521739130434,
            ),
            (
                "GSE165870 bulk support",
                "score_GSE165870_bulk_support",
                0.06,
                0.04695652173913043,
            ),
            (
                "Network pharmacology support",
                "score_network_pharmacology_hit",
                0.20,
                0.21739130434782608,
            ),
        ],
        columns=[
            "evidence_component",
            "score_column",
            "pre_removal_weight",
            "final_normalized_weight",
        ],
    )
    derivation["derivation"] = (
        "The archived docking weight (0.08) was removed; the remaining "
        "experimental total (0.72) and network-pharmacology weight (0.20) "
        "were divided by 0.92 without changing their relative proportions."
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    scenario_out.to_csv(
        args.output_dir / "top30_weight_sensitivity_scenarios.csv", index=False
    )
    selected_out.to_csv(
        args.output_dir / "top30_weight_sensitivity_selected_long.csv", index=False
    )
    gene_out.to_csv(
        args.output_dir / "top30_weight_sensitivity_gene_frequency.csv", index=False
    )
    summary_out.to_csv(
        args.output_dir / "top30_weight_sensitivity_summary.csv", index=False
    )
    derivation.to_csv(args.output_dir / "top30_weight_derivation.csv", index=False)

    primary = scenario_out.loc[
        scenario_out["scenario_id"].eq("primary_no_docking")
    ]
    if len(primary) != 1:
        raise ValueError("Exactly one primary_no_docking scenario is required.")
    if primary.iloc[0]["top30_genes"].split(", ") != formal_order:
        raise ValueError(
            "Primary scenario order does not reproduce the frozen formal Top30."
        )

    print(
        json.dumps(
            {
                "eligible_candidates": int(len(eligible)),
                "scenarios": int(len(scenario_out)),
                "primary_top30_reproduced": True,
                "minimum_overlap_count": int(
                    scenario_out["top30_overlap_count"].min()
                ),
                "median_overlap_count": float(
                    scenario_out["top30_overlap_count"].median()
                ),
                "output_dir": str(args.output_dir),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
