from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


QC = Path(r"E:\AA新投稿思路\BMC_Pharmacology_and_Toxicology_投稿_20260830\90_内部QC_勿提交")
OUTPUT = QC / "simulated_editor_sensitivity_20260902_manifest_sha256_v1.csv"
if OUTPUT.exists():
    raise FileExistsError(f"Refusing to overwrite {OUTPUT}")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


patterns = [
    "simulated_editor_sensitivity_20260902_v2_*",
    "simulated_editor_baseline_only_wgcna_20260902_v1_*",
    "run_simulated_editor_sensitivity_20260902_v1.py",
    "run_baseline_only_wgcna_20260902_v1.R",
    "summarize_baseline_only_wgcna_mapping_20260902_v1.py",
    "simulated_editor_sensitivity_20260902_execution_record.md",
]
files = sorted({path for pattern in patterns for path in QC.glob(pattern) if path.is_file()})
rows = [
    {
        "file": path.name,
        "bytes": path.stat().st_size,
        "sha256": digest(path),
    }
    for path in files
]
pd.DataFrame(rows).to_csv(OUTPUT, index=False)
print(f"Wrote {len(rows)} records to {OUTPUT}")
