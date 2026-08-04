#!/usr/bin/env python3
"""Validate the packaged formal scTenifoldKnk design and candidate outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


CORE10 = ["CDK6", "CA2", "PARP1", "KIT", "SYK", "GSK3B", "HIF1A", "TOP2A", "TERT", "CD38"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    args = parser.parse_args()
    data_dir = args.data_dir
    manifest_path = data_dir / "run_manifest_210.csv"
    endpoints_path = data_dir / "descriptive_endpoints_all_runs.csv"
    diff_dir = data_dir / "current210_candidate_diffregulation_20260725"
    for path in [manifest_path, endpoints_path, diff_dir]:
        if not path.exists():
            raise FileNotFoundError(path)

    manifest = pd.read_csv(manifest_path)
    endpoints = pd.read_csv(endpoints_path)
    if len(manifest) != 210 or len(endpoints) != 210:
        raise ValueError("Expected exactly 210 manifest and endpoint rows")
    if manifest["run_id"].duplicated().any() or endpoints["run_id"].duplicated().any():
        raise ValueError("Duplicate run_id detected")
    if set(manifest["run_id"]) != set(endpoints["run_id"]):
        raise ValueError("Manifest and endpoint run identifiers differ")
    if set(manifest["common_seed"].astype(int)) != {20260724}:
        raise ValueError("Formal common seed must equal 20260724 in all 210 runs")

    roles = manifest["run_role"].value_counts().to_dict()
    if roles.get("candidate") != 10 or roles.get("matched_control") != 200:
        raise ValueError(f"Unexpected role counts: {roles}")
    candidates = manifest.loc[manifest["run_role"].eq("candidate"), "gKO"].tolist()
    if candidates != CORE10:
        raise ValueError(f"Candidate order mismatch: {candidates}")
    controls = manifest.loc[manifest["run_role"].eq("matched_control"), "gKO"]
    if controls.duplicated().any():
        raise ValueError("Matched controls are not disjoint")

    files = sorted(diff_dir.glob("*_candidate_*_diffRegulation.csv"))
    if len(files) != 10:
        raise ValueError("Expected 10 candidate differential-response files")
    observed = []
    for path in files:
        frame = pd.read_csv(path)
        if list(frame.columns) != ["gene", "distance", "Z", "FC", "p.value", "p.adj"]:
            raise ValueError(f"Unexpected columns in {path.name}")
        if len(frame) != 3000:
            raise ValueError(f"Expected 3000 genes in {path.name}")
        observed.append(path.name.split("_candidate_", 1)[1].split("_diffRegulation", 1)[0])
    if observed != CORE10:
        raise ValueError(f"Differential-response candidate order mismatch: {observed}")

    print("PASS: formal 10 + 200 = 210 scTenifoldKnk design verified")


if __name__ == "__main__":
    main()
