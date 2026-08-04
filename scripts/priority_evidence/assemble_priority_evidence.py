#!/usr/bin/env python3
"""Assemble the prespecified 30-candidate evidence annotation without re-ranking.

This entry point joins author-staged evidence tables to a pre-specified candidate
order. Docking, molecular-dynamics, and perturbation outputs are rejected as
inputs and are never used for selection or ranking.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


FORBIDDEN_COLUMNS = re.compile(r"docking|vina|molecular.?dynamics|perturb|rank|score_total", re.I)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def gene_key(row: dict[str, str]) -> str:
    for key in ("GeneSymbol", "gene", "gene_symbol", "symbol"):
        if row.get(key, "").strip():
            return row[key].strip().upper()
    raise ValueError("A GeneSymbol/gene/gene_symbol/symbol column is required")


def index_rows(rows: list[dict[str, str]], label: str) -> dict[str, dict[str, str]]:
    index: dict[str, dict[str, str]] = {}
    for row in rows:
        key = gene_key(row)
        if key in index:
            raise ValueError(f"Duplicate gene {key!r} in {label}")
        index[key] = row
    return index


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-order", type=Path, required=True,
                        help="TSV containing the prespecified 30-gene order")
    parser.add_argument("--evidence-table", type=Path, required=True,
                        help="TSV of non-docking evidence annotations")
    parser.add_argument("--wgcna-table", type=Path, required=True,
                        help="TSV of strict42 module-level annotations")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--qc-json", type=Path, required=True)
    args = parser.parse_args()

    for output in (args.output, args.qc_json):
        if output.exists():
            raise SystemExit(f"Refusing to overwrite existing output: {output}")

    order_rows = read_tsv(args.candidate_order)
    evidence_rows = read_tsv(args.evidence_table)
    wgcna_rows = read_tsv(args.wgcna_table)
    order = [gene_key(row) for row in order_rows]
    if len(order) != 30 or len(set(order)) != 30:
        raise ValueError("The prespecified candidate order must contain 30 unique genes")

    for label, rows in (("evidence", evidence_rows), ("WGCNA", wgcna_rows)):
        fields = set().union(*(row.keys() for row in rows)) if rows else set()
        rejected = sorted(field for field in fields if FORBIDDEN_COLUMNS.search(field))
        if rejected:
            raise ValueError(f"Forbidden post-prioritization columns in {label}: {rejected}")

    evidence = index_rows(evidence_rows, "evidence table")
    wgcna = index_rows(wgcna_rows, "WGCNA table")
    missing_evidence = [gene for gene in order if gene not in evidence]
    missing_wgcna = [gene for gene in order if gene not in wgcna]
    if missing_evidence or missing_wgcna:
        raise ValueError(
            f"Missing rows: evidence={missing_evidence}; WGCNA={missing_wgcna}"
        )

    evidence_fields = [f for f in evidence_rows[0] if f not in {"GeneSymbol", "gene", "gene_symbol", "symbol"}]
    wgcna_fields = [f for f in wgcna_rows[0] if f not in {"GeneSymbol", "gene", "gene_symbol", "symbol"}]
    output_fields = ["fixed_order", "GeneSymbol"] + [f"evidence__{f}" for f in evidence_fields] + [f"wgcna__{f}" for f in wgcna_fields]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()
        for rank, gene in enumerate(order, start=1):
            row = {"fixed_order": rank, "GeneSymbol": gene}
            row.update({f"evidence__{f}": evidence[gene].get(f, "") for f in evidence_fields})
            row.update({f"wgcna__{f}": wgcna[gene].get(f, "") for f in wgcna_fields})
            writer.writerow(row)

    qc = {
        "candidate_count": len(order),
        "candidate_order_preserved": True,
        "ranking_performed": False,
        "docking_used": False,
        "molecular_dynamics_used": False,
        "computational_perturbation_used": False,
        "missing_evidence": missing_evidence,
        "missing_wgcna": missing_wgcna,
    }
    args.qc_json.write_text(json.dumps(qc, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

