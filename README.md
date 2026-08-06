# Donor-aware recovery and perturbation workflow for aplastic anemia

Private public-release candidate containing analysis code, selected derived
source tables, and a frozen publication-replay test suite for the BMC Genomics
manuscript. The repository covers prescription
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
- `reference_outputs/`: frozen submission figures used only as regression-test references.
- `FIGURE_SOURCE_MAP.tsv`: Figure 1-9 and supplementary source mapping.
- `REPRODUCIBILITY_MATRIX.tsv`: publication replay versus scientific recomputation.
- `PUBLICATION_ASSET_CHECKSUMS.tsv`: frozen v5 submission-asset checksums.
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

## Publication replay

Reproducibility is reported in two layers and the distinction is mandatory:

1. **Publication replay** rebuilds the final Figure 1-9 layouts, Figure S8-S10,
   and the eight-page Figure S1-S8 package from repository-distributed derived
   tables and approved publication intermediates. It then compares the results
   with frozen BMC submission references.
2. **Scientific recomputation** reruns analytical models from public,
   author-staged, or licensed inputs. Some workflows require separately
   obtained GEO files, reviewed design metadata, official model assets,
   prepared structural inputs, or licensed database exports.

The publication replay is available on Windows with Python 3.10+ and R 4.5.x:

Install the dedicated replay stack from
`environment/publication_replay_python_20260806.txt`; do not substitute the
Geneformer GPU environment for this rendering-only command.

```powershell
py -3 -B scripts/publication_figures/reproduce_all_publication_figures.py `
  --output-root C:\AA_replay_20260806 `
  --rscript "C:\Program Files\R\R-4.5.2\bin\Rscript.exe"
```

The output root must not already exist. Figures 1, 2, 4, 5, 6, 8, 9, S9 and
S10 require exact SHA-256 equality. Figures 3, 7 and S8 use documented pixel
tolerances because supported Matplotlib/R rasterizer versions and the frozen
submission re-encoding change antialiasing or PNG metadata without changing
source values or panel geometry.

To verify a separately held frozen BMC submission package, including the
corrected 12-rule supplementary workbook, run:

```powershell
py -3 -B scripts/publication_tables/verify_submission_assets.py `
  C:\path\to\BMC_Genomics_Submission_Ready_20260806_v5_Pending_GitHub `
  --inspect-workbooks
```

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

Licensed database exports, prepared docking structures, model weights,
official third-party source snapshots, private participant design metadata,
and raw/topology/trajectory MD files are excluded. The journal supplementary
workbooks are verified as separate publication assets rather than duplicated
into GitHub; this avoids treating archived database rows in the submission
workbook as a license to redistribute the underlying source export. Included
derived tables may retain public-study pseudonymous subject labels and
single-cell barcodes for provenance; they contain no names or contact details.
Public GEO accessions and exact staging requirements are listed in
`DATA_AND_LICENSES.md` and `config/input_paths.example.tsv`.

## Figure reproducibility

The complete publication appearance of Figure 1-9 and Figure S1-S10 is
replayable or regression-verifiable. This does **not** mean that every analysis
is self-contained from raw data. Figures 1, 2, 4, 5 and 6 and several legacy
supplementary panels still require author-staged or licensed upstream inputs
for scientific recomputation. Figure 7 is rerendered from the complete derived
summary of five 100 ns systems (500 ns total), while raw trajectories remain
excluded. Figure 9/S9/S10 is rerendered from frozen model outputs; repeating
model inference requires official Geneformer assets. Exact statuses are in
`REPRODUCIBILITY_MATRIX.tsv` and `FIGURE_SOURCE_MAP.tsv`.

The finalized Figure 9/S9/S10 presentation pipeline uses the v6 renderer,
followed by the v7 RGB conversion. These versions preserve black titles and a
white background; v7 changes color mode only and does not recompute results.

## License and citation

Repository-authored code is BSD-3-Clause. Third-party datasets, databases,
Geneformer code/model assets, SPSS, Cytoscape, GROMACS, AutoDock Vina and
R/Python packages retain their own terms. The repository URL is recorded in
`CITATION.cff`; a Zenodo DOI must be added only after one is actually issued.

Software name: `aa-donor-recovery-perturbation`; supported operating systems:
Windows for the byte-stable publication replay and Windows/Linux for analytical
scripts where dependencies permit; languages: Python 3.10+ and R 4.5.x; license:
BSD-3-Clause for repository-authored code. External software, data and model
restrictions are documented in `DATA_AND_LICENSES.md`.
