from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"E:\AA新投稿思路")
QC = ROOT / "BMC_Pharmacology_and_Toxicology_投稿_20260830" / "90_内部QC_勿提交"
REPO = QC / "github_sync_worktree"
HSPC = (
    ROOT
    / "02_项目总归档_拍板版"
    / "04_前置方法与算法"
    / "方案前置方法"
    / "02_GEO_WGCNA_原始数据与设计"
    / "04_GSE247531_CD34_HSPC分析结果"
)
LEGACY_TABLES = (
    ROOT
    / "02_项目总归档_拍板版"
    / "02_正式文稿与汇报"
    / "02_所有结果汇总审阅"
    / "68_Pharmaceuticals投稿前故事线整理_20260630"
    / "05_tables"
    / "legacy_20260610_current_method_tables"
)

SCORE_INPUT = LEGACY_TABLES / "81_Tier1_core_targets_no_docking_formal_scores.csv"
META_INPUT = HSPC / "92_CD34_WGCNA_strict_donor_aware_20260715_collapsed_sample_metadata.tsv"
EXPR_INPUT = HSPC / "92_CD34_WGCNA_strict_donor_aware_20260715_collapsed_log1pCP10K_by_subject_timepoint.csv"
DONOR_SCORE_INPUT = REPO / "derived_data" / "geneformer" / "AA_candidate_set_directional_recovery_extension_MVP_20260802_v1_donor_scores.csv"
MATCH_DIAG_INPUT = (
    ROOT
    / "02_项目总归档_拍板版"
    / "01_拍板研究方案"
    / "BMC_Genomics_PRISM_Select_MVP_matched_null_diagnostics_20260801_v2.csv"
)
MD_ROOT = REPO / "derived_data" / "molecular_dynamics"

PREFIX = "simulated_editor_sensitivity_20260902_v2"
OUTPUTS = {
    "ranking_long": QC / f"{PREFIX}_outcome_blind_ranking_long.csv",
    "ranking_summary": QC / f"{PREFIX}_outcome_blind_summary.csv",
    "trajectory": QC / f"{PREFIX}_six_month_trajectory.csv",
    "trajectory_summary": QC / f"{PREFIX}_trajectory_summary.csv",
    "null_distribution": QC / f"{PREFIX}_trajectory_null_distribution.csv.gz",
    "null_summary": QC / f"{PREFIX}_trajectory_null_summary.csv",
    "md_inventory": QC / f"{PREFIX}_md_qc_inventory.csv",
    "md_coverage": QC / f"{PREFIX}_md_system_coverage.csv",
    "report": QC / f"{PREFIX}_report.md",
    "runtime": QC / f"{PREFIX}_runtime.json",
}

FORMAL_WEIGHTS = {
    "score_WGCNA_module_strength": 0.2347826087,
    "score_module_pathology_importance": 0.1408695652,
    "score_CD34_expression": 0.1408695652,
    "score_CD34_marker_context": 0.1095652174,
    "score_AA_direction_recovery": 0.1095652174,
    "score_GSE165870_bulk_support": 0.0469565217,
    "score_network_pharmacology_hit": 0.2173913043,
}

SCENARIOS = {
    "formal_all_components": [],
    "drop_explicit_recovery_only": ["score_AA_direction_recovery"],
    "drop_recovery_and_pathology_trait": [
        "score_AA_direction_recovery",
        "score_module_pathology_importance",
    ],
    "drop_all_module_and_recovery_features": [
        "score_WGCNA_module_strength",
        "score_module_pathology_importance",
        "score_AA_direction_recovery",
    ],
}

CORE10 = ["CDK6", "PARP1", "KIT", "SYK", "HIF1A", "TOP2A", "CA2", "CD38", "TERT", "GSK3B"]
FOCUS5 = ["TOP2A", "GSK3B", "KIT", "HIF1A", "SYK"]
N_PERM = 2000
SEED = 20260902


def refuse_overwrite() -> None:
    existing = [str(path) for path in OUTPUTS.values() if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing outputs: " + "; ".join(existing))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def exact_two_sided_sign_p(n_positive: int, n_total: int) -> float:
    tail_n = min(n_positive, n_total - n_positive)
    tail = sum(math.comb(n_total, k) for k in range(tail_n + 1))
    return min(1.0, 2.0 * tail / (2**n_total))


def markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(column) for column in frame.columns]
    header = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join(["---"] * len(columns)) + " |"
    rows = []
    for row in frame.itertuples(index=False, name=None):
        values = []
        for value in row:
            if isinstance(value, float):
                values.append(f"{value:.6g}")
            else:
                values.append(str(value))
        rows.append("| " + " | ".join(values) + " |")
    return "\n".join([header, divider, *rows])


def bootstrap_median(values: np.ndarray, seed: int, n_boot: int = 2000) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    samples = rng.choice(values, size=(n_boot, len(values)), replace=True)
    medians = np.median(samples, axis=1)
    return tuple(np.quantile(medians, [0.025, 0.975]).tolist())


def rank_with_weights(df: pd.DataFrame, dropped: list[str]) -> pd.DataFrame:
    weights = {key: value for key, value in FORMAL_WEIGHTS.items() if key not in dropped}
    scale = sum(weights.values())
    if scale <= 0:
        raise ValueError("No nonzero ranking weights remain")
    score = sum(df[column].astype(float) * weight for column, weight in weights.items()) / scale * 100.0
    result = df[["GeneSymbol", "rank_no_docking_formal", "no_docking_formal_score"]].copy()
    result["scenario_score"] = score
    result = result.sort_values(["scenario_score", "GeneSymbol"], ascending=[False, True]).reset_index(drop=True)
    result["scenario_rank"] = np.arange(1, len(result) + 1)
    return result


def outcome_blind_analysis() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    started = time.perf_counter()
    df = pd.read_csv(SCORE_INPUT)
    for column in FORMAL_WEIGHTS:
        df[column] = pd.to_numeric(df[column], errors="raise")

    formal = rank_with_weights(df, [])
    merged_formal = formal.merge(df[["GeneSymbol", "no_docking_formal_score"]], on="GeneSymbol", suffixes=("", "_source"))
    max_abs_diff = float(np.max(np.abs(merged_formal["scenario_score"] - merged_formal["no_docking_formal_score_source"])))
    if max_abs_diff > 1e-6:
        raise ValueError(f"Formal score reproduction failed; max absolute difference={max_abs_diff}")

    formal_rank = formal.set_index("GeneSymbol")["scenario_rank"]
    formal_top30 = set(formal.head(30)["GeneSymbol"])
    formal_top10 = set(formal.head(10)["GeneSymbol"])
    long_parts: list[pd.DataFrame] = []
    summaries: list[dict] = []
    for scenario, dropped in SCENARIOS.items():
        ranked = rank_with_weights(df, dropped)
        ranked.insert(0, "scenario", scenario)
        ranked["formal_rank_reproduced"] = ranked["GeneSymbol"].map(formal_rank)
        ranked["rank_shift_vs_formal"] = ranked["scenario_rank"] - ranked["formal_rank_reproduced"]
        ranked["in_scenario_top30"] = ranked["scenario_rank"] <= 30
        ranked["is_core10"] = ranked["GeneSymbol"].isin(CORE10)
        ranked["is_focus5"] = ranked["GeneSymbol"].isin(FOCUS5)
        ranked["dropped_components"] = ";".join(dropped) if dropped else "none"
        long_parts.append(ranked)

        current_top30 = set(ranked.head(30)["GeneSymbol"])
        current_top10 = set(ranked.head(10)["GeneSymbol"])
        focus_ranks = ranked.set_index("GeneSymbol")["scenario_rank"].to_dict()
        summaries.append(
            {
                "scenario": scenario,
                "dropped_components": ";".join(dropped) if dropped else "none",
                "top30_overlap_n": len(current_top30 & formal_top30),
                "top30_jaccard": len(current_top30 & formal_top30) / len(current_top30 | formal_top30),
                "top10_overlap_n": len(current_top10 & formal_top10),
                "spearman_rank_vs_formal": float(
                    np.corrcoef(
                        ranked["scenario_rank"].to_numpy(dtype=float),
                        ranked["formal_rank_reproduced"].to_numpy(dtype=float),
                    )[0, 1]
                ),
                "core10_retained_top30_n": sum(focus_ranks[g] <= 30 for g in CORE10),
                **{f"rank_{g}": int(focus_ranks[g]) for g in FOCUS5},
            }
        )
    ranking_long = pd.concat(long_parts, ignore_index=True)
    ranking_summary = pd.DataFrame(summaries)
    runtime = {"seconds": time.perf_counter() - started, "formal_score_max_abs_diff": max_abs_diff, "input_rows": len(df)}
    return ranking_long, ranking_summary, runtime


def prepare_recovery_folds() -> tuple[pd.DataFrame, list[dict]]:
    meta = pd.read_csv(META_INPUT, sep="\t")
    expr = pd.read_csv(EXPR_INPUT)
    sample_ids = meta["subject_timepoint_id"].tolist()
    matrix = expr.set_index("gene")[sample_ids].T.astype(float)
    matrix.index.name = "subject_timepoint_id"
    meta = meta.set_index("subject_timepoint_id").loc[sample_ids]
    hd_ids = meta.index[meta["disease"].eq("HD")].tolist()
    baseline_meta = meta[(meta["disease"].eq("SAA")) & (meta["timepoint"].eq("baseline"))]
    match_diag = pd.read_csv(MATCH_DIAG_INPUT)
    primary_match = match_diag[match_diag["k_controls"].eq(20)]
    excluded_predictor_genes = set(CORE10)
    for value in primary_match["selected_control_genes"].fillna(""):
        excluded_predictor_genes.update(gene for gene in str(value).split(";") if gene)
    records: list[dict] = []
    folds: list[dict] = []

    for n_top_genes in (1000, 3000, 5000):
        for subject in sorted(baseline_meta["subject"].unique()):
            subject_rows = meta[(meta["subject"].eq(subject)) & (meta["disease"].eq("SAA"))]
            baseline_ids = subject_rows.index[subject_rows["timepoint"].eq("baseline")].tolist()
            follow_ids = subject_rows.index[subject_rows["timepoint"].isin(["3M", "6M"])].tolist()
            if len(baseline_ids) != 1 or not follow_ids:
                continue
            held_baseline = baseline_ids[0]
            train_baselines = baseline_meta.index[baseline_meta["subject"].ne(subject)].tolist()
            anchor_ids = hd_ids + train_baselines
            anchor = matrix.loc[anchor_ids]
            gene_var = anchor.var(axis=0, ddof=1).replace([np.inf, -np.inf], np.nan).dropna()
            gene_var = gene_var[(gene_var > 1e-10) & (~gene_var.index.isin(excluded_predictor_genes))]
            selected_genes = gene_var.nlargest(min(n_top_genes, len(gene_var))).index
            anchor_selected = anchor[selected_genes].to_numpy(dtype=float)
            mu = anchor_selected.mean(axis=0)
            sd = anchor_selected.std(axis=0, ddof=1)
            usable = np.isfinite(sd) & (sd > 1e-10)
            selected_genes = selected_genes[usable]
            mu = mu[usable]
            sd = sd[usable]
            z = (matrix.loc[:, selected_genes] - mu) / sd
            hd_centroid = z.loc[hd_ids].mean(axis=0).to_numpy()
            base_centroid = z.loc[train_baselines].mean(axis=0).to_numpy()
            axis = hd_centroid - base_centroid
            denominator = float(np.dot(axis, axis))
            baseline_vec = z.loc[held_baseline].to_numpy()

            def score(sample_id: str) -> float:
                return float(np.dot(z.loc[sample_id].to_numpy() - base_centroid, axis) / denominator)

            for follow_id in follow_ids:
                timepoint = str(meta.loc[follow_id, "timepoint"])
                baseline_score = score(held_baseline)
                followup_score = score(follow_id)
                records.append(
                    {
                        "gene_set_size": n_top_genes,
                        "subject": subject,
                        "followup_timepoint": timepoint,
                        "baseline_id": held_baseline,
                        "followup_id": follow_id,
                        "baseline_recovery_score": baseline_score,
                        "followup_recovery_score": followup_score,
                        "paired_recovery_shift": followup_score - baseline_score,
                        "n_genes_in_fold": len(selected_genes),
                    }
                )
                if n_top_genes == 5000:
                    folds.append(
                        {
                            "subject": subject,
                            "followup_timepoint": timepoint,
                            "axis": axis,
                            "denominator": denominator,
                            "difference": z.loc[follow_id].to_numpy() - baseline_vec,
                            "anchor_z": z.loc[anchor_ids].to_numpy(),
                            "n_true_hd": len(hd_ids),
                        }
                    )
    return pd.DataFrame(records), folds


def trajectory_analysis() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    started = time.perf_counter()
    trajectory, all_folds = prepare_recovery_folds()
    trajectory["time_month"] = trajectory["followup_timepoint"].map({"3M": 3, "6M": 6}).astype(int)
    trajectory["is_latest_available"] = trajectory["time_month"] == trajectory.groupby(["gene_set_size", "subject"])["time_month"].transform("max")

    source_scores = pd.read_csv(DONOR_SCORE_INPUT)
    expected = source_scores[source_scores["test_id"].eq("primary_all10_equal")][
        ["subject", "followup_timepoint", "independent_recovery_shift"]
    ].drop_duplicates()
    check = trajectory[(trajectory["gene_set_size"].eq(5000)) & trajectory["is_latest_available"]].merge(
        expected, on=["subject", "followup_timepoint"], how="left"
    )
    max_abs_diff = float(np.max(np.abs(check["paired_recovery_shift"] - check["independent_recovery_shift"])))
    if max_abs_diff > 1e-10:
        raise ValueError(f"Recovery-axis reproduction failed; max absolute difference={max_abs_diff}")

    summaries: list[dict] = []
    for label, block in {
        "latest_available_5000": trajectory[(trajectory["gene_set_size"].eq(5000)) & trajectory["is_latest_available"]],
        "six_month_only_5000": trajectory[(trajectory["gene_set_size"].eq(5000)) & trajectory["followup_timepoint"].eq("6M")],
        "six_month_only_3000": trajectory[(trajectory["gene_set_size"].eq(3000)) & trajectory["followup_timepoint"].eq("6M")],
        "six_month_only_1000": trajectory[(trajectory["gene_set_size"].eq(1000)) & trajectory["followup_timepoint"].eq("6M")],
    }.items():
        values = block["paired_recovery_shift"].to_numpy(dtype=float)
        n_pos = int(np.sum(values > 0))
        ci_low, ci_high = bootstrap_median(values, SEED + int(block["gene_set_size"].iloc[0]))
        summaries.append(
            {
                "analysis": label,
                "n_donors": len(values),
                "n_positive": n_pos,
                "positive_fraction": n_pos / len(values),
                "median_shift": float(np.median(values)),
                "median_bootstrap_95ci_low": ci_low,
                "median_bootstrap_95ci_high": ci_high,
                "exact_two_sided_sign_p": exact_two_sided_sign_p(n_pos, len(values)),
                "minimum_shift": float(np.min(values)),
                "maximum_shift": float(np.max(values)),
            }
        )
    trajectory_summary = pd.DataFrame(summaries)

    latest_subjects = set(
        trajectory[(trajectory["gene_set_size"].eq(5000)) & trajectory["is_latest_available"]]["subject"]
    )
    folds = []
    for fold in all_folds:
        if fold["subject"] not in latest_subjects:
            continue
        subject_latest = trajectory[(trajectory["gene_set_size"].eq(5000)) & trajectory["subject"].eq(fold["subject"]) & trajectory["is_latest_available"]]
        if len(subject_latest) != 1 or fold["followup_timepoint"] != subject_latest.iloc[0]["followup_timepoint"]:
            continue
        folds.append(fold)
    if len(folds) != 17:
        raise ValueError(f"Expected 17 latest-available folds, found {len(folds)}")

    rng_axis = np.random.default_rng(SEED + 101)
    rng_label = np.random.default_rng(SEED + 202)
    observed_shifts = np.array([np.dot(f["difference"], f["axis"]) / f["denominator"] for f in folds])
    observed_median = float(np.median(observed_shifts))
    observed_positive = int(np.sum(observed_shifts > 0))
    null_rows: list[dict] = []
    for iteration in range(1, N_PERM + 1):
        random_axis_shifts = []
        label_shifts = []
        for fold in folds:
            perm_axis = fold["axis"][rng_axis.permutation(len(fold["axis"]))]
            random_axis_shifts.append(float(np.dot(fold["difference"], perm_axis) / np.dot(perm_axis, perm_axis)))

            anchor_z = fold["anchor_z"]
            n_hd = fold["n_true_hd"]
            pseudo_hd_idx = rng_label.choice(anchor_z.shape[0], size=n_hd, replace=False)
            mask = np.ones(anchor_z.shape[0], dtype=bool)
            mask[pseudo_hd_idx] = False
            null_axis = anchor_z[pseudo_hd_idx].mean(axis=0) - anchor_z[mask].mean(axis=0)
            null_denom = float(np.dot(null_axis, null_axis))
            label_shifts.append(float(np.dot(fold["difference"], null_axis) / null_denom))

        for null_type, values in (
            ("gene_weight_permuted_axis", np.asarray(random_axis_shifts)),
            ("anchor_label_permutation", np.asarray(label_shifts)),
        ):
            null_rows.append(
                {
                    "null_type": null_type,
                    "iteration": iteration,
                    "median_shift": float(np.median(values)),
                    "n_positive": int(np.sum(values > 0)),
                }
            )
    null_distribution = pd.DataFrame(null_rows)
    null_summaries: list[dict] = []
    for null_type, block in null_distribution.groupby("null_type"):
        null_summaries.append(
            {
                "null_type": null_type,
                "iterations": len(block),
                "observed_median_shift": observed_median,
                "null_median_shift_mean": float(block["median_shift"].mean()),
                "null_median_shift_q025": float(block["median_shift"].quantile(0.025)),
                "null_median_shift_q975": float(block["median_shift"].quantile(0.975)),
                "one_sided_empirical_p_median": float((1 + np.sum(block["median_shift"] >= observed_median)) / (len(block) + 1)),
                "observed_n_positive": observed_positive,
                "null_n_positive_q025": float(block["n_positive"].quantile(0.025)),
                "null_n_positive_q975": float(block["n_positive"].quantile(0.975)),
                "one_sided_empirical_p_n_positive": float((1 + np.sum(block["n_positive"] >= observed_positive)) / (len(block) + 1)),
            }
        )
    runtime = {
        "seconds": time.perf_counter() - started,
        "axis_reproduction_max_abs_diff": max_abs_diff,
        "permutations_per_null": N_PERM,
        "latest_donors": len(folds),
    }
    return trajectory, trajectory_summary, null_distribution, pd.DataFrame(null_summaries), runtime


def count_table_rows_and_columns(path: Path) -> tuple[int | None, int | None, str]:
    suffixes = path.suffixes
    is_gzip = suffixes and suffixes[-1].lower() == ".gz"
    data_suffix = suffixes[-2].lower() if is_gzip and len(suffixes) >= 2 else path.suffix.lower()
    if data_suffix not in {".csv", ".tsv"}:
        return None, None, ""
    opener = gzip.open if is_gzip else open
    delimiter = "," if data_suffix == ".csv" else "\t"
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        try:
            header = next(reader)
        except StopIteration:
            return 0, 0, ""
        rows = sum(1 for _ in reader)
    return rows, len(header), ";".join(header)


def md_inventory() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    started = time.perf_counter()
    rows: list[dict] = []
    for path in sorted(MD_ROOT.rglob("*")):
        if not path.is_file():
            continue
        n_rows, n_columns, columns = count_table_rows_and_columns(path)
        rows.append(
            {
                "relative_path": path.relative_to(REPO).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "table_rows": n_rows,
                "table_columns": n_columns,
                "columns": columns,
            }
        )
    inventory = pd.DataFrame(rows)

    summary_path = MD_ROOT / "current_figure8" / "Figure8_final20ns_quantitative_summary.tsv"
    systems = pd.read_csv(summary_path, sep="\t")[["case_id", "display"]].copy()
    systems["core_100ns_timeseries_present"] = True
    systems["ca_rmsf_present"] = True
    systems["final20ns_summary_present"] = True
    systems["thermodynamic_qc_temperature_pressure_density_present"] = systems["case_id"].eq("TOP2A_sesamin")
    systems["protein_sasa_qc_present"] = systems["case_id"].eq("TOP2A_sesamin")
    systems["paired_apo_comparator_present"] = systems["case_id"].eq("HIF1A_ARNT_butin_exploratory")
    systems["extended_interface_pocket_qc_present"] = systems["case_id"].eq("HIF1A_ARNT_butin_exploratory")
    apo = pd.DataFrame(
        [
            {
                "case_id": "HIF1A_ARNT_apo",
                "display": "HIF1A–ARNT apo",
                "core_100ns_timeseries_present": True,
                "ca_rmsf_present": True,
                "final20ns_summary_present": False,
                "thermodynamic_qc_temperature_pressure_density_present": False,
                "protein_sasa_qc_present": True,
                "paired_apo_comparator_present": True,
                "extended_interface_pocket_qc_present": True,
            }
        ]
    )
    coverage = pd.concat([systems, apo], ignore_index=True)
    coverage["complete_raw_trajectory_archived_in_repository"] = False
    coverage["replicate_trajectories_present"] = False
    coverage["interpretive_limit"] = "single-trajectory descriptive support; no binding free energy or inferential claim"
    runtime = {"seconds": time.perf_counter() - started, "file_count": len(inventory), "system_count": len(coverage)}
    return inventory, coverage, runtime


def main() -> None:
    refuse_overwrite()
    wall_start = time.perf_counter()
    clock_start = time.strftime("%Y-%m-%dT%H:%M:%S%z")

    ranking_long, ranking_summary, ranking_runtime = outcome_blind_analysis()
    trajectory, trajectory_summary, null_distribution, null_summary, trajectory_runtime = trajectory_analysis()
    inventory, coverage, md_runtime = md_inventory()

    ranking_long.to_csv(OUTPUTS["ranking_long"], index=False)
    ranking_summary.to_csv(OUTPUTS["ranking_summary"], index=False)
    trajectory.to_csv(OUTPUTS["trajectory"], index=False)
    trajectory_summary.to_csv(OUTPUTS["trajectory_summary"], index=False)
    null_distribution.to_csv(OUTPUTS["null_distribution"], index=False, compression="gzip")
    null_summary.to_csv(OUTPUTS["null_summary"], index=False)
    inventory.to_csv(OUTPUTS["md_inventory"], index=False)
    coverage.to_csv(OUTPUTS["md_coverage"], index=False)

    sixm = trajectory_summary.set_index("analysis").loc["six_month_only_5000"]
    strict = ranking_summary.set_index("scenario").loc["drop_recovery_and_pathology_trait"]
    ultra = ranking_summary.set_index("scenario").loc["drop_all_module_and_recovery_features"]
    report = f"""# Simulated-editor sensitivity analyses (internal QC; not submission-ready by itself)

Generated: {time.strftime('%Y-%m-%d %H:%M:%S %z')}

## Scope and safeguards

- No formal manuscript, figure, supplement, GitHub branch, Zenodo record, or raw molecular-dynamics trajectory was changed.
- The formal 126-gene score table was reproduced before any component was removed (maximum absolute score difference {ranking_runtime['formal_score_max_abs_diff']:.3g}).
- The 5,000-gene donor recovery scores were independently reconstructed from the collapsed subject-by-timepoint expression matrix and matched the archived donor scores (maximum absolute difference {trajectory_runtime['axis_reproduction_max_abs_diff']:.3g}).
- All sensitivity outputs are newly written under `90_内部QC_勿提交`; scripts refuse to overwrite existing outputs.

## 1. Outcome-blind ranking sensitivity

The nested analyses are deliberately labeled by what was removed rather than presented as a single post hoc replacement ranking.

- Removing the explicit longitudinal AA-direction/recovery component retained {int(ranking_summary.set_index('scenario').loc['drop_explicit_recovery_only','top30_overlap_n'])}/30 formal top-30 genes.
- Removing both explicit recovery and the module-pathology trait component retained {int(strict['top30_overlap_n'])}/30 formal top-30 genes (Jaccard {strict['top30_jaccard']:.3f}; rank Spearman {strict['spearman_rank_vs_formal']:.3f}).
- The more conservative analysis that removed WGCNA module strength as well as pathology and recovery retained {int(ultra['top30_overlap_n'])}/30 formal top-30 genes (Jaccard {ultra['top30_jaccard']:.3f}; rank Spearman {ultra['spearman_rank_vs_formal']:.3f}).
- Focus-five ranks under the recovery-plus-pathology-blind analysis: TOP2A {int(strict['rank_TOP2A'])}, GSK3B {int(strict['rank_GSK3B'])}, KIT {int(strict['rank_KIT'])}, HIF1A {int(strict['rank_HIF1A'])}, SYK {int(strict['rank_SYK'])}.

Interpretation boundary: this is a robustness analysis of the prespecified score components. It does not create independent biological validation, and the strictest scenario still retains cell-intrinsic expression/marker-context features derived from the same public cohort.

## 2. Six-month-only healthy-directed trajectory

- The 6-month-only primary 5,000-gene analysis included {int(sixm['n_donors'])} donors; {int(sixm['n_positive'])}/{int(sixm['n_donors'])} shifted in the healthy-directed direction.
- Median displacement was {sixm['median_shift']:.6f} (donor-bootstrap 95% interval {sixm['median_bootstrap_95ci_low']:.6f} to {sixm['median_bootstrap_95ci_high']:.6f}); exact two-sided sign-test P = {sixm['exact_two_sided_sign_p']:.8g}.
- The 1,000- and 3,000-gene six-month sensitivity specifications are reported in the companion summary table.

## 3. Random-axis and anchor-label null tests

Each null used {N_PERM} deterministic Monte Carlo replicates across the 17 latest-available donor folds. Gene-weight permutation preserves each fold's observed axis norm but scrambles its gene alignment; anchor-label permutation preserves the 4-versus-18 group sizes within each fold.

{markdown_table(null_summary)}

These empirical tests assess whether the observed aggregate projection is unusually aligned with the prespecified healthy-versus-AA direction. They are not substitutes for an external cohort.

## 4. Six-system molecular-dynamics QC inventory

- Five ligand-containing systems have archived per-frame core metrics, C-alpha RMSF, and final-20-ns summaries.
- HIF1A-ARNT apo is archived as the sixth comparator with extended interface/pocket metrics.
- TOP2A has temperature, pressure, density, SASA, coordinate, and contact-occupancy QC tables.
- HIF1A-ARNT apo/butin has paired interface, chain, pocket, SASA, hydrogen-bond, and RMSF tables.
- GSK3B, KIT, and SYK do not have system-specific thermodynamic (temperature/pressure/density) QC tables in the current public derived-data tree.
- No system has replicate trajectories or complete raw trajectories archived in this repository; therefore the MD evidence remains descriptive and should not be phrased as binding-affinity or inferential validation.

## 5. Baseline-only WGCNA status

The exact collapsed input contains 23 one-profile-per-participant baseline/reference profiles (19 SAA baseline and 4 healthy donors), so a baseline-only network is computationally feasible today. However, four healthy donors are too few for strong between-group module-trait inference. A separate no-overwrite R run is used only to quantify module reproducibility and feasibility; it must not replace the 42-profile donor-aware sensitivity without explicit interpretation of its lower power.
"""
    OUTPUTS["report"].write_text(report, encoding="utf-8")

    runtime = {
        "clock_started": clock_start,
        "clock_finished": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "total_seconds": time.perf_counter() - wall_start,
        "ranking": ranking_runtime,
        "trajectory_and_nulls": trajectory_runtime,
        "md_inventory": md_runtime,
        "inputs": {str(path): {"bytes": path.stat().st_size, "sha256": sha256(path)} for path in [SCORE_INPUT, META_INPUT, EXPR_INPUT, DONOR_SCORE_INPUT, MATCH_DIAG_INPUT]},
        "outputs": {key: str(path) for key, path in OUTPUTS.items()},
    }
    OUTPUTS["runtime"].write_text(json.dumps(runtime, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(runtime, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
