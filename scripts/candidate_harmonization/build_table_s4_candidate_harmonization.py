#!/usr/bin/env python3
"""Build Supplementary Table S4 from archived, separately sourced inputs.

TCMSP supplies the retained herb-compound inventory. Archived
SwissTargetPrediction output supplies compound-target records. The two sources
are never conflated. Composite target labels are retained for provenance and
excluded from exact single-symbol matching.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path

import pandas as pd


EXPECTED_DATABASE_COUNTS = {
    "TTD": 2,
    "OMIM": 10,
    "DISGENET": 4,
    "GENE CARD": 1528,
}


def norm(value: object) -> str:
    return str(value if value is not None else "").replace("\ufeff", "").strip()


def header_key(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", norm(value).lower())


def find_column(columns, candidates) -> str:
    keyed = {header_key(column): column for column in columns}
    for candidate in candidates:
        key = header_key(candidate)
        if key in keyed:
            return keyed[key]
    for candidate in candidates:
        key = header_key(candidate)
        for observed, column in keyed.items():
            if key in observed or observed in key:
                return column
    raise ValueError(f"Required column not found: {candidates}")


def herb_name(value: object) -> str:
    text = norm(value)
    key = text.lower()
    if re.search(r"eclipta|prostrata|hanliancao|mohanlian|mo han lian|墨旱莲", key):
        return "Ecliptae Herba"
    if re.search(r"cuscuta|chinensis|tusizi|tu si zi|菟丝子", key):
        return "Cuscutae Semen"
    if re.search(r"ligustr|lucidum|nüzhenzi|nvzhenzi|nu zhen zi|女贞子", key):
        return "Ligustri Lucidi Fructus"
    raise ValueError(f"Unrecognized herb label: {text}")


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest().upper()


def read_first_sheet(path: Path) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name=0, dtype=str).dropna(how="all")


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")


def yn(value: bool) -> str:
    return "Yes" if value else "No"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory-xlsx", type=Path, required=True)
    parser.add_argument("--target-raw-csv", type=Path, required=True)
    parser.add_argument("--target-dedup-csv", type=Path, required=True)
    parser.add_argument("--disease-raw-csv", type=Path, required=True)
    parser.add_argument("--disease-dedup-csv", type=Path, required=True)
    parser.add_argument("--intersection-csv", type=Path, required=True)
    parser.add_argument("--output-xlsx", type=Path, required=True)
    args = parser.parse_args()

    inputs = [
        args.inventory_xlsx,
        args.target_raw_csv,
        args.target_dedup_csv,
        args.disease_raw_csv,
        args.disease_dedup_csv,
        args.intersection_csv,
    ]
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.output_xlsx.exists():
        raise FileExistsError(f"Refusing to overwrite: {args.output_xlsx}")
    args.output_xlsx.parent.mkdir(parents=True, exist_ok=True)

    inventory = read_first_sheet(args.inventory_xlsx)
    inv_herb = find_column(inventory.columns, ["Herb name", "Herb", "Drug"])
    inv_mol = find_column(inventory.columns, ["MOL ID", "MOLID", "molecule_id"])
    inv_compound = find_column(
        inventory.columns, ["Molecule name", "Compound name", "Active component"]
    )
    inventory = inventory.loc[inventory[inv_mol].map(norm).ne("")].copy()
    inventory["herb_source_label"] = inventory[inv_herb].map(norm)
    inventory["herb_standard_name"] = inventory[inv_herb].map(herb_name)
    inventory["mol_id"] = inventory[inv_mol].map(norm)
    inventory["compound_name_as_recorded"] = inventory[inv_compound].map(norm)

    target = read_csv(args.target_raw_csv)
    target_herb = find_column(target.columns, ["drug", "herb"])
    target_mol = find_column(target.columns, ["MOLID", "MOL ID"])
    target_compound = find_column(target.columns, ["moleculename", "compound name"])
    target_entity = find_column(target.columns, ["GeneName", "target entity"])
    target = target.loc[target[target_entity].map(norm).ne("")].copy()
    target["herb_source_label"] = target[target_herb].map(norm)
    target["herb_standard_name"] = target[target_herb].map(herb_name)
    target["mol_id"] = target[target_mol].map(norm)
    target["compound_name_as_recorded"] = target[target_compound].map(norm)
    target["source_target_entity"] = target[target_entity].map(norm)
    target["entity_type"] = target["source_target_entity"].map(
        lambda x: "composite_target_entity" if re.search(r"\s", x) else "single_gene_symbol_form"
    )
    relation_cols = [
        "herb_source_label",
        "mol_id",
        "compound_name_as_recorded",
        "source_target_entity",
    ]
    target["exact_relation_occurrence_n"] = (
        target.groupby(relation_cols, sort=False)["source_target_entity"].transform("size")
    )
    target["first_exact_occurrence"] = target.groupby(
        relation_cols, sort=False
    ).cumcount().eq(0).map(yn)
    target.insert(0, "source_row_id", range(1, len(target) + 1))

    target_dedup = read_csv(args.target_dedup_csv)
    source_entities = [norm(value) for value in target_dedup.iloc[:, 0] if norm(value)]
    disease_raw = read_csv(args.disease_raw_csv)
    disease_sets = {
        norm(column).upper(): {norm(value) for value in disease_raw[column] if norm(value)}
        for column in disease_raw.columns
    }
    disease_dedup = read_csv(args.disease_dedup_csv)
    aa_genes = [norm(value) for value in disease_dedup.iloc[:, 0] if norm(value)]
    intersection = read_csv(args.intersection_csv)
    candidates = [norm(value) for value in intersection.iloc[:, 0] if norm(value)]

    relation_unique = target[relation_cols].drop_duplicates()
    entity_set = set(target["source_target_entity"])
    single_entities = sorted(
        value for value in entity_set if not re.search(r"\s", value)
    )
    composite_entities = sorted(value for value in entity_set if re.search(r"\s", value))
    aa_set = set(aa_genes)
    candidate_set = set(candidates)
    exact_intersection = set(single_entities) & aa_set

    checks = {
        "inventory_rows": len(inventory) == 34,
        "unique_inventory_mol": inventory["mol_id"].nunique() == 29,
        "target_rows": len(target) == 2357,
        "distinct_relations": len(relation_unique) == 2344,
        "target_entities": len(entity_set) == 481,
        "source_entity_crosscheck": len(source_entities) == 481
        and set(source_entities) == entity_set,
        "single_entities": len(single_entities) == 472,
        "composite_entities": len(composite_entities) == 9,
        "aa_genes": len(aa_genes) == 1529 and len(aa_set) == 1529,
        "candidate_genes": len(candidates) == 126 and len(candidate_set) == 126,
        "exact_intersection": exact_intersection == candidate_set,
        "database_counts": all(
            len(disease_sets.get(name, set())) == count
            for name, count in EXPECTED_DATABASE_COUNTS.items()
        ),
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError("Frozen source checks failed: " + ", ".join(failed))

    entity_rows = []
    for entity in source_entities:
        rows = target.loc[target["source_target_entity"].eq(entity)]
        entity_rows.append(
            {
                "source_target_entity": entity,
                "entity_type": (
                    "composite_target_entity"
                    if re.search(r"\s", entity)
                    else "single_gene_symbol_form"
                ),
                "retained_relationship_rows_n": len(rows),
                "distinct_exact_relationships_n": len(rows[relation_cols].drop_duplicates()),
                "herb_n": rows["herb_standard_name"].nunique(),
                "compound_n": rows["mol_id"].nunique(),
                "candidate_pool_exact_match": yn(entity in candidate_set and not re.search(r"\s", entity)),
            }
        )
    entity_frame = pd.DataFrame(entity_rows)

    def memberships(gene: str) -> dict[str, object]:
        flags = {
            f"in_{name.replace(' ', '_')}_retained_list": yn(gene in values)
            for name, values in disease_sets.items()
        }
        flags["database_membership_n"] = sum(
            gene in values for values in disease_sets.values()
        )
        return flags

    aa_frame = pd.DataFrame(
        [{"gene_symbol": gene, **memberships(gene), "candidate_pool_exact_match": yn(gene in candidate_set)}
         for gene in aa_genes]
    )
    candidate_frame = pd.DataFrame(
        [
            {
                "candidate_index": index,
                "gene_symbol": gene,
                **memberships(gene),
                "retained_relationship_rows_n": int(
                    target["source_target_entity"].eq(gene).sum()
                ),
                "distinct_exact_relationships_n": len(
                    target.loc[target["source_target_entity"].eq(gene), relation_cols]
                    .drop_duplicates()
                ),
                "matching_rule": "Exact standardized single-gene-symbol match",
            }
            for index, gene in enumerate(candidates, 1)
        ]
    )

    inventory_key_counts = target.groupby(
        ["herb_standard_name", "mol_id"], sort=False
    ).size()
    inventory_out = pd.DataFrame(
        {
            "inventory_row_id": range(1, len(inventory) + 1),
            "herb_source_label": inventory["herb_source_label"],
            "herb_standard_name": inventory["herb_standard_name"],
            "mol_id": inventory["mol_id"],
            "compound_name_as_recorded": inventory["compound_name_as_recorded"],
        }
    )
    inventory_out["retained_relationship_rows_n"] = [
        int(inventory_key_counts.get((herb, mol), 0))
        for herb, mol in zip(
            inventory_out["herb_standard_name"], inventory_out["mol_id"]
        )
    ]
    inventory_out["has_archived_target_relations"] = (
        inventory_out["retained_relationship_rows_n"].gt(0).map(yn)
    )

    summary = pd.DataFrame(
        [
            ("herb_compound_records", 34),
            ("unique_MOL_identifiers", 29),
            ("retained_relationship_rows", 2357),
            ("distinct_exact_relationships", 2344),
            ("source_target_entities", 481),
            ("single_gene_symbol_forms", 472),
            ("composite_target_entities", 9),
            ("AA_resource_genes", 1529),
            ("candidate_target_pool", 126),
        ],
        columns=["metric", "value"],
    )
    readme = pd.DataFrame(
        [
            ("Scope", "TCMSP compound inventory, archived SwissTargetPrediction records, AA database genes, and the exact 126-gene candidate pool."),
            ("Boundary", "The exact intersection defines a candidate space and does not establish binding, disease causality, or therapeutic efficacy."),
            ("Matching", "Only single standardized gene-symbol entities enter exact matching; nine composite labels remain visible but are not decomposed."),
        ],
        columns=["item", "description"],
    )
    source_index = pd.DataFrame(
        [
            {
                "source_id": f"S4-SRC-{index:02d}",
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "scientific_role": role,
            }
            for index, (path, role) in enumerate(
                zip(
                    inputs,
                    [
                        "TCMSP compound inventory",
                        "archived SwissTargetPrediction target records",
                        "target-entity cross-check",
                        "four AA database gene lists",
                        "exact AA gene-symbol union",
                        "exact candidate-pool list",
                    ],
                ),
                1,
            )
        ]
    )

    with pd.ExcelWriter(args.output_xlsx, engine="xlsxwriter") as writer:
        sheets = {
            "README": readme,
            "S4_summary": summary,
            "S4_herb_compounds": inventory_out,
            "S4_target_relations": target[
                [
                    "source_row_id",
                    "herb_source_label",
                    "herb_standard_name",
                    "mol_id",
                    "compound_name_as_recorded",
                    "source_target_entity",
                    "entity_type",
                    "exact_relation_occurrence_n",
                    "first_exact_occurrence",
                ]
            ],
            "S4_target_entities": entity_frame,
            "S4_AA_resource_genes": aa_frame,
            "S4_candidate_pool": candidate_frame,
            "Source_index": source_index,
        }
        for sheet_name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=sheet_name, index=False)
            sheet = writer.sheets[sheet_name]
            sheet.freeze_panes(1, 0)
            sheet.autofilter(0, 0, max(len(frame), 1), max(len(frame.columns) - 1, 0))
            for column_index, column in enumerate(frame.columns):
                width = min(
                    48,
                    max(len(str(column)) + 2, frame[column].astype(str).map(len).max() + 2),
                )
                sheet.set_column(column_index, column_index, width)

    print(f"PASS: wrote {args.output_xlsx}")


if __name__ == "__main__":
    main()
