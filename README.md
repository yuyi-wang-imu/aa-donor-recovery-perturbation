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

This repository is currently a **private, pre-release candidate** at
<https://github.com/yuyi-wang-imu/aa-donor-recovery-perturbation>. It is not a
GitHub release and has no Zenodo DOI. The historical repository and historical
DOI are not modified or withdrawn by this work and must not be reused as the
citation for this repository.

## Repository map

- `scripts/run_workflow.py`: lists, checks, and dispatches canonical entries.
- `WORKFLOW_ORDER.tsv`: intended execution order and workflow scope.
- `config/`: parameters and external-input staging template.
- `environment/`: recorded software versions and install specifications.
- `derived_data/`: selected figure source tables; no raw trajectories or model weights.
- `FIGURE_SOURCE_MAP.tsv`: Figure 1-9 and supplementary mapping, including gaps.
- `DATA_AND_LICENSES.md`: redistribution boundaries and third-party assets.
- `CITATION.cff`: citation metadata without a guessed DOI.
- `MANIFEST.tsv`: canonical path, size and SHA-256 inventory.

## Quick structural checks

Python 3.10 or newer is required. On Linux/macOS, run:

```bash
python3 scripts/run_workflow.py --list
python3 scripts/run_workflow.py --check
python3 scripts/validate_repository.py
```

On Windows, where `python` may resolve to a legacy interpreter, use:

```powershell
py -3 -B scripts/run_workflow.py --list
py -3 -B scripts/run_workflow.py --check
py -3 -B scripts/validate_repository.py
```

These checks verify repository structure; they do not download licensed data,
Geneformer weights, or large molecular-dynamics trajectories and do not rerun
the scientific analyses.

## Environment

The supported non-Geneformer package ranges are in
`environment/requirements.txt`; R dependencies are in
`environment/r_packages.tsv`. The host snapshot in
`environment/python_packages.tsv` is provenance, not a dependency lock. See
`environment/README.md` for the distinction. The recorded GPU Geneformer run
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
private participant design metadata, and raw/topology/trajectory MD files are
excluded. Included derived tables may retain public-study pseudonymous subject
labels and single-cell barcodes for provenance; they contain no names or contact
details. Public GEO accessions and exact staging requirements are listed in
`DATA_AND_LICENSES.md` and `config/input_paths.example.tsv`.

## Figure reproducibility

Figures 2, 3, 7, 8, 9 and the mapped supplementary Geneformer panels have
packaged source tables and/or canonical renderers. Figures 1, 4, 5, 6 and parts
of the legacy supplementary package depend on accepted images or external
inputs that are not fully redistributed. These are explicitly labelled in
`FIGURE_SOURCE_MAP.tsv`; absence of an exact renderer must not be described as
full from-raw reproducibility.

The finalized Figure 9/S9/S10 presentation pipeline uses the v6 renderer,
followed by the v7 RGB conversion. These versions preserve black titles and a
white background; v7 changes color mode only and does not recompute results.

## License and citation

Repository-authored code is BSD-3-Clause. Third-party datasets, databases,
Geneformer code/model assets, SPSS, Cytoscape, GROMACS, AutoDock Vina and
R/Python packages retain their own terms. The repository URL is recorded in
`CITATION.cff`; a Zenodo DOI must be added only after one is actually issued.
