from __future__ import annotations

from pathlib import Path

import pandas as pd


QC = Path(r"E:\AA新投稿思路\BMC_Pharmacology_and_Toxicology_投稿_20260830\90_内部QC_勿提交")
PREFIX = "simulated_editor_baseline_only_wgcna_20260902_v1"
CANDIDATE_INPUT = QC / f"{PREFIX}_candidate10_mapping.csv"
BEST_PAIR_INPUT = QC / f"{PREFIX}_original42_overlap_best_pairs.csv"
MAPPING_OUTPUT = QC / f"{PREFIX}_candidate10_best_overlap_mapping.csv"
ADDENDUM_OUTPUT = QC / f"{PREFIX}_interpretive_addendum.md"

existing = [str(path) for path in (MAPPING_OUTPUT, ADDENDUM_OUTPUT) if path.exists()]
if existing:
    raise FileExistsError("Refusing to overwrite existing outputs: " + "; ".join(existing))

candidates = pd.read_csv(CANDIDATE_INPUT)
best = pd.read_csv(BEST_PAIR_INPUT).rename(
    columns={
        "original42_module": "best_overlap_original42_module",
        "jaccard": "best_overlap_jaccard",
        "overlap_n": "best_overlap_n",
    }
)
keep = [
    "baseline_module",
    "best_overlap_original42_module",
    "best_overlap_jaccard",
    "best_overlap_n",
]
result = candidates.merge(best[keep], on="baseline_module", how="left")
result["best_overlap_recovers_candidate_original_module"] = (
    result["best_overlap_original42_module"] == result["original42_module"]
)
result.to_csv(MAPPING_OUTPUT, index=False)

covered = result[result["baseline_module"].notna() & result["original42_module"].notna()]
n_match = int(covered["best_overlap_recovers_candidate_original_module"].sum())
addendum = f"""# Interpretive addendum: candidate module correspondence

WGCNA color labels are nominal, so exact color-name agreement is not a valid module-preservation measure by itself. Mapping each baseline-only module to its best-overlap original 42-profile module showed that {n_match}/{len(covered)} expression-covered frozen candidates recovered their original module by best Jaccard correspondence. CDK6, PARP1, KIT, SYK, HIF1A, TOP2A, and CA2 matched; CD38 did not. GSK3B and TERT were not in the baseline-only top-5,000 MAD network and therefore cannot be evaluated in this mapping.

This supports partial module preservation after removing repeated longitudinal profiles, but the overall adjusted Rand index (0.3044) indicates only moderate global concordance. The 23-profile sensitivity remains underpowered because it contains only four healthy donors.
"""
ADDENDUM_OUTPUT.write_text(addendum, encoding="utf-8")
print(result.to_string(index=False))
