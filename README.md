# Donor-aware recovery and perturbation workflow for aplastic anemia

Public-release candidate containing analysis code and selected derived source
tables for the BMC Genomics manuscript. The repository covers prescription
mining, exact candidate construction, transcriptomic and co-expression
analyses, bone-marrow single-cell projection, evidence integration, docking,
five standardized 100 ns molecular-dynamics simulations (500 ns total),
scTenifoldKnk perturbation, and donor-aware Geneformer deletion and
overexpression analyses.

The manuscript emphasizes the stronger donor-aware bidirectional support for
TOP2A and GSK3B, treats KIT cautiously, and keeps complete weak/null calibration
results in Supplementary Table S10. Model outputs are computational evidence,
not experimental target validation.

## Release status

This directory is an **unpublished candidate**. It is not a GitHub release and
has no new DOI. The historical repository and historical Zenodo DOI are not
modified or withdrawn by this work and must not be reused as the citation for
this candidate. See `UPLOAD_READINESS_AUDIT_20260804_v1.md` before uploading.

## Repository map

- `scripts/run_workflow.py`: lists, checks, and dispatches canonical entries.
- `WORKFLOW_ORDER.tsv`: intended execution order and workflow scope.
- `config/`: parameters and external-input staging template.
- `environment/`: recorded software versions and install specifications.
- `derived_data/`: selected figure source tables; no raw trajectories or model weights.
- `FIGURE_SOURCE_MAP.tsv`: Figure 1-9 and supplementary mapping, including gaps.
- `DATA_AND_LICENSES.md`: redistribution boundaries and third-party assets.
- `PROPOSED_UPLOAD_MANIFEST_20260804_v1.tsv`: upload/include decision by path group.

## Quick structural checks

```bash
python scripts/run_workflow.py --list
python scripts/run_workflow.py --check
python scripts/validate_repository.py
```

These checks verify repository structure; they do not download licensed data,
Geneformer weights, or large molecular-dynamics trajectories and do not rerun
the scientific analyses.

## Environment

The original non-Geneformer workflows use `environment/requirements.txt` and
the R packages in `environment/r_packages.tsv`. The recorded GPU Geneformer run
used Python 3.10.20, NumPy 1.26.4, pandas 2.2.3, SciPy 1.12.0,
PyTorch 2.4.1+cu118, Transformers 4.46.0, and CUDA. A reproducible starting
specification is provided in
`environment/geneformer_gpu_environment_20260804_v1.yml`.

GPU execution accelerates model inference but does not constitute a new
algorithm by itself. The publishable methodological contribution is the
donor-aware recovery-axis design, bidirectional single-gene perturbation,
donor-bootstrap inference, state stratification, and matched-background
calibration. CPU execution may be technically possible but was not the recorded
runtime and has not been clean-room validated for this release.

## Geneformer assets and inputs

1. Obtain the public GSE247531 input from NCBI GEO and an author-reviewed design
   table with the columns documented in `config/input_paths.example.tsv`.
2. Run the balanced-input preparation script with a new output directory:

```bash
Rscript scripts/geneformer/70_bmc_geneformer_prepare_balanced_cd34_public_20260804_v1.R \
  /path/to/GSE247531_CD34counts.Rdata.gz \
  /path/to/GSE247531_CD34_design_table.tsv \
  /path/to/balanced_cd34_inputs
```

3. Obtain the official Geneformer source dictionaries and Geneformer-V2-104M
   model separately. Set `GENEFORMER_SOURCE_DIR` and `GENEFORMER_MODEL_DIR`.
   The model weights are intentionally excluded because the recorded file is
   approximately 418 MB and is maintained by the upstream project.
4. Run the entries in `WORKFLOW_ORDER.tsv`. Each analysis must write to a new
   output directory; never point a public command at the archived project tree.

The recorded model SHA-256 was
`fff5cba29ddd8792991fa77b4872246fbe548a178cebda3775cdc72b67780e7f`.

## Inputs not redistributed

Author-curated prescription records, licensed database exports, prepared
docking structures, model weights, official third-party source snapshots,
participant design metadata pending privacy review, and raw/topology/trajectory
MD files are excluded. Public GEO accessions and exact staging requirements are
listed in `DATA_AND_LICENSES.md` and `config/input_paths.example.tsv`.

## Figure reproducibility

Figures 2, 3, 7, 8, 9 and the mapped supplementary Geneformer panels have
packaged source tables and/or canonical renderers. Figures 1, 4, 5, 6 and parts
of the legacy supplementary package depend on accepted images or external
inputs that are not fully redistributed. These are explicitly labelled in
`FIGURE_SOURCE_MAP.tsv`; absence of an exact renderer must not be described as
full from-raw reproducibility.

The final publication choice between the visually reviewed Figure 9/S9/S10 v5
set and later v6/v7 presentation derivatives remains an author decision.

## License and citation

Repository-authored code is BSD-3-Clause. Third-party datasets, databases,
Geneformer code/model assets, SPSS, Cytoscape, GROMACS, AutoDock Vina and R/Python
packages retain their own terms. Do not add a new GitHub URL or Zenodo DOI to
`CITATION.cff` until the author creates the new repository and release record.
