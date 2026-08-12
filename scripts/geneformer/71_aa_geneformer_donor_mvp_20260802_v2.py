from __future__ import annotations

import importlib.util
from pathlib import Path


V1 = Path(__file__).with_name("71_aa_geneformer_donor_mvp_20260802_v1.py")
spec = importlib.util.spec_from_file_location("aa_geneformer_mvp_v1", V1)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Unable to load frozen MVP implementation: {V1}")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def find_preferred_source_asset(pattern: str) -> Path:
    preferred = module.SOURCE_DIR / "geneformer" / pattern
    if preferred.exists():
        return preferred
    hits = [path for path in sorted(module.SOURCE_DIR.rglob(pattern)) if "/build/lib/" not in path.as_posix()]
    if len(hits) != 1:
        raise RuntimeError(f"Expected one non-build source asset for {pattern}; observed {len(hits)}: {hits}")
    return hits[0]


module.find_one = find_preferred_source_asset

if __name__ == "__main__":
    raise SystemExit(module.main())
