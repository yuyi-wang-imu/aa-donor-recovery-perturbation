#!/usr/bin/env python3
"""Run the pre-specified 6 x 23 AutoDock Vina matrix without filtering.

The design TSV is the sole source of receptor, ligand, and pocket information.
All runs and scores are retained. No affinity threshold or target ranking is
implemented. The default mode performs preflight only; use --execute to run.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from pathlib import Path


REQUIRED = {
    "run_id", "target", "ligand", "receptor_pdbqt", "ligand_pdbqt",
    "center_x", "center_y", "center_z",
}


def parse_mode1(log_path: Path) -> str:
    if not log_path.exists():
        return ""
    pattern = re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+")
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = pattern.match(line)
        if match:
            return match.group(1)
    return ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--vina", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-runs", type=int, default=138)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()

    if args.output_dir.exists():
        raise SystemExit(f"Refusing to overwrite existing output directory: {args.output_dir}")

    with args.design.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not REQUIRED.issubset(rows[0]):
        raise ValueError(f"Design columns must include {sorted(REQUIRED)}")
    if len(rows) != args.expected_runs:
        raise ValueError(f"Expected {args.expected_runs} runs, found {len(rows)}")
    if len({row["run_id"] for row in rows}) != len(rows):
        raise ValueError("run_id values must be unique")

    planned = []
    for row in rows:
        receptor = Path(row["receptor_pdbqt"])
        ligand = Path(row["ligand_pdbqt"])
        planned.append({
            **row,
            "receptor_exists": str(receptor.is_file()),
            "ligand_exists": str(ligand.is_file()),
        })
    if args.execute and not args.vina.is_file():
        raise FileNotFoundError(args.vina)
    missing = [row["run_id"] for row in planned if row["receptor_exists"] != "True" or row["ligand_exists"] != "True"]
    if args.execute and missing:
        raise FileNotFoundError(f"Missing receptor/ligand inputs for runs: {missing}")

    args.output_dir.mkdir(parents=True)
    results = []
    for row in planned:
        run_dir = args.output_dir / row["run_id"]
        command = [
            str(args.vina),
            "--receptor", row["receptor_pdbqt"],
            "--ligand", row["ligand_pdbqt"],
            "--center_x", row["center_x"],
            "--center_y", row["center_y"],
            "--center_z", row["center_z"],
            "--size_x", "15", "--size_y", "15", "--size_z", "15",
            "--exhaustiveness", "8", "--num_modes", "9",
            "--out", str(run_dir / "pose.pdbqt"),
            "--log", str(run_dir / "vina.log"),
        ]
        status = "planned"
        returncode = ""
        if args.execute:
            run_dir.mkdir()
            proc = subprocess.run(command, capture_output=True, text=True, check=False)
            (run_dir / "stdout_stderr.txt").write_text(proc.stdout + "\n" + proc.stderr, encoding="utf-8")
            returncode = str(proc.returncode)
            status = "completed" if proc.returncode == 0 else "failed"
        results.append({
            "run_id": row["run_id"], "target": row["target"], "ligand": row["ligand"],
            "status": status, "returncode": returncode,
            "mode1_affinity_kcal_mol": parse_mode1(run_dir / "vina.log") if args.execute else "",
            "command": json.dumps(command, ensure_ascii=False),
        })

    with (args.output_dir / "all_runs.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    (args.output_dir / "run_parameters.json").write_text(json.dumps({
        "expected_runs": args.expected_runs,
        "box_angstrom": [15, 15, 15],
        "exhaustiveness": 8,
        "num_modes": 9,
        "seed": "not prespecified",
        "filtering_or_ranking": False,
        "executed": args.execute,
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

