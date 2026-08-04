#!/usr/bin/env python3
"""Generate a PBC-aware GROMACS trajectory-analysis plan or execute it.

The design TSV supplies explicit group selections for each system. The script
produces trajectory metrics only and does not label a complex as stable.
Default behavior is a dry run; use --execute after validating group names.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path


REQUIRED = {
    "system_id", "tpr", "xtc", "index", "center_group", "output_group",
    "backbone_group", "calpha_group", "protein_group", "ligand_group",
}


def run(command: list[str], selections: str, cwd: Path) -> dict:
    proc = subprocess.run(command, input=selections, text=True, capture_output=True, cwd=cwd, check=False)
    return {"command": command, "selections": selections.strip().splitlines(), "returncode": proc.returncode,
            "stdout": proc.stdout, "stderr": proc.stderr}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()

    if args.output_dir.exists():
        raise SystemExit(f"Refusing to overwrite existing output directory: {args.output_dir}")
    with args.design.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not REQUIRED.issubset(rows[0]):
        raise ValueError(f"Design columns must include {sorted(REQUIRED)}")
    if len({row["system_id"] for row in rows}) != len(rows):
        raise ValueError("system_id values must be unique")

    args.output_dir.mkdir(parents=True)
    plan = []
    for row in rows:
        system_dir = args.output_dir / row["system_id"]
        tpr, xtc, ndx = Path(row["tpr"]), Path(row["xtc"]), Path(row["index"])
        inputs_exist = all(path.is_file() for path in (tpr, xtc, ndx))
        pbc_xtc = system_dir / "trajectory_pbc_centered.xtc"
        commands = [
            ([args.gmx, "trjconv", "-s", str(tpr), "-f", str(xtc), "-n", str(ndx), "-o", str(pbc_xtc), "-pbc", "mol", "-center", "-ur", "compact"], f"{row['center_group']}\n{row['output_group']}\n"),
            ([args.gmx, "rms", "-s", str(tpr), "-f", str(pbc_xtc), "-n", str(ndx), "-o", str(system_dir / "backbone_rmsd.xvg"), "-tu", "ns"], f"{row['backbone_group']}\n{row['backbone_group']}\n"),
            ([args.gmx, "rmsf", "-s", str(tpr), "-f", str(pbc_xtc), "-n", str(ndx), "-o", str(system_dir / "calpha_rmsf.xvg"), "-res"], f"{row['calpha_group']}\n"),
            ([args.gmx, "gyrate", "-s", str(tpr), "-f", str(pbc_xtc), "-n", str(ndx), "-o", str(system_dir / "protein_rg.xvg")], f"{row['protein_group']}\n"),
            ([args.gmx, "mindist", "-s", str(tpr), "-f", str(pbc_xtc), "-n", str(ndx), "-od", str(system_dir / "protein_ligand_min_distance.xvg")], f"{row['protein_group']}\n{row['ligand_group']}\n"),
        ]
        entry = {"system_id": row["system_id"], "inputs_exist": inputs_exist, "executed": args.execute, "commands": []}
        if args.execute:
            if not inputs_exist:
                raise FileNotFoundError(f"Missing input for {row['system_id']}")
            system_dir.mkdir()
            for command, selections in commands:
                result = run(command, selections, system_dir)
                entry["commands"].append(result)
                if result["returncode"] != 0:
                    raise RuntimeError(f"GROMACS command failed for {row['system_id']}: {result}")
        else:
            entry["commands"] = [{"command": command, "selections": selections.strip().splitlines()} for command, selections in commands]
        plan.append(entry)

    (args.output_dir / "trajectory_analysis_plan.json").write_text(json.dumps({
        "pbc_strategy": "whole-molecule reconstruction followed by centering in a compact unit cell",
        "time_axis": "read from the trajectory and reported in ns where supported",
        "interpretation_note": "trajectory metrics are descriptive and do not establish experimental binding",
        "systems": plan,
    }, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
