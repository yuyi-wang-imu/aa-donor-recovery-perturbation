# Donor-aware longitudinal transcriptomics and bidirectional Geneformer perturbation in aplastic anemia

This repository contains reproducible analysis code, configuration files,
selected derived source tables, and publication-regression tests supporting the
manuscript *Donor-Aware Longitudinal Transcriptomics and Bidirectional
Geneformer Perturbation Prioritize Recovery-Aligned Genes in Aplastic Anemia*.
The workflow separates traceable candidate provenance from recovery-aligned
model evidence. It integrates prescription-derived candidate construction,
CD34+ hematopoietic transcriptomics, participant-aware co-expression analysis,
bone-marrow single-cell context, donor-aware longitudinal recovery-axis
analysis, paired in silico deletion and overexpression with Geneformer,
scTenifoldKnk, docking, and molecular dynamics.

## Evidence hierarchy and current manuscript alignment

The current manuscript follows this prespecified evidence order:

1. Mining 390 aplastic-anemia prescription records identified a recurrent
   *Ecliptae Herba*-*Cuscutae Semen*-*Ligustri Lucidi Fructus* pattern.
2. Exact matching of 472 standardized compound-associated genes to 1,529
   aplastic-anemia-associated genes defined 126 starting candidates.
3. CD34+ expression, participant-aware WGCNA, bone-marrow single-cell context,
   independent expression context, and prespecified multi-source criteria
   prioritized 30 candidates and selected ten genes before model perturbation.
4. All 17 donors with pretreatment and latest-follow-up profiles moved toward
   the healthy reference on the primary 5,000-gene recovery axis (median
   displacement, 0.362; exact two-sided sign-test P = 1.53 x 10^-5).
5. Geneformer reproduced the healthy-directed change in 14 of 15
   baseline-to-six-month donor pairs. Five candidates--TOP2A, KIT,
   GSK3B, HIF1A, and SYK--showed the prespecified deletion-away and
   overexpression-toward-health pattern; TOP2A and GSK3B received the strongest
   donor-level and expression-matched support.

scTenifoldKnk, docking, and molecular dynamics are downstream,
hypothesis-generating analyses. They did not contribute to candidate selection
or ranking.
Prescription recurrence is not evidence of efficacy; Geneformer output is
model-based prioritization rather than experimental gene editing; docking
scores are not biochemical affinities; and single 100 ns trajectories do not
establish direct target engagement, inhibition, or therapeutic activity.

## Version and publication-asset boundary

Release `v0.1.0` is the immutable software snapshot archived at Zenodo and
corresponds to commit
`b850d7e6a6141962d3a2d38021d0780e3c5907a1`. The DOI and tag are retained for
traceability. Maintenance commits on the default branch, including the
2026-08-20 manuscript-structure alignment, do not alter the archived ZIP, tag,
or DOI.

The current Human Genomics manuscript uses eight main figures: CD34+/WGCNA;
bone-marrow single-cell context; candidate annotation and sensitivity; the
three-herb/compound/candidate network; donor recovery plus Geneformer;
scTenifoldKnk; representative docking hypotheses; and molecular dynamics.
Supplementary Figures run from S1 to S16; the original prescription-mining and
candidate-construction workflow is retained in Figure S11. The current
structural follow-up reports five ligand-containing 100 ns trajectories plus a
matched 100 ns HIF1A-ARNT apo reference (600 ns total).
The archived replay tables and reference images distributed in `v0.1.0` retain
their original Figure 1-9 and Figure S1-S10 numbering. Files such as
`FIGURE_SOURCE_MAP.tsv`, `REPRODUCIBILITY_MATRIX.tsv`, and
`PUBLICATION_ASSET_CHECKSUMS.tsv` are therefore archive-specific regression
records, not a claim that the current manuscript still uses the legacy figure
numbering.

## Repository status

The source code and supporting files are maintained at
<https://github.com/yuyi-wang-imu/aa-donor-recovery-perturbation>. Release
`v0.1.0` is archived at <https://doi.org/10.5281/zenodo.21837457>.

## Repository map

- `scripts/run_workflow.py`: lists, checks, and dispatches canonical entries.
- `WORKFLOW_ORDER.tsv`: intended execution order and workflow scope.
- `config/`: parameters and external-input staging template.
- `environment/`: recorded software versions and install specifications.
- `derived_data/`: selected figure source tables; no raw trajectories or model weights.
- `reference_outputs/`: publication figures used as regression-test references.
- `reference_outputs/current_manuscript/`: the eight current main figures and
  graphical abstract, with exact checksums.
- `CURRENT_MANUSCRIPT_FIGURE_SOURCE_MAP.tsv`: current Figure 1-8 source and
  access boundaries.
- `CURRENT_MANUSCRIPT_REPRODUCIBILITY_MATRIX.tsv`: current publication replay
  versus scientific recomputation.
- `CURRENT_MANUSCRIPT_ASSET_CHECKSUMS.tsv`: exact 19-file Human Genomics
  submission-package inventory used by the current verifier.
- `CURRENT_MANUSCRIPT_SUPPLEMENTARY_FIGURE_MAP.tsv`: page-level Figure S1-S16
  content and verification scope.
- `FIGURE_SOURCE_MAP.tsv`: archive-specific Figure 1-9 and supplementary
  source mapping for the archived `v0.1.0` replay.
- `REPRODUCIBILITY_MATRIX.tsv`: publication replay versus scientific recomputation.
- `PUBLICATION_ASSET_CHECKSUMS.tsv`: checksums for the nine main figures and
  six Additional files in the archived `v0.1.0` submission-asset snapshot.
- `DATA_AND_LICENSES.md`: redistribution boundaries and third-party assets.
- `CITATION.cff`: citation metadata including the issued Zenodo DOI.
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

## Current-manuscript verification and publication replay

Reproducibility is reported in four layers and the distinction is mandatory:

1. **Current submission verification** checks the separately held flat
   19-file Human Genomics package against exact byte sizes and SHA-256 values,
   verifies Figure 1-8 dimensions and resolution metadata, verifies the
   920 x 300-pixel graphical abstract and 16-page supplementary PDF, and can
   inspect core workbook invariants. The package is not copied into GitHub.
2. **Current Figure 8 replay** rerenders the five-candidate 100 ns comparison
   from the repository-distributed time series, C-alpha RMSF and final-20-ns
   summaries. Raw trajectories, topologies and checkpoints are excluded.
3. **Archived publication replay** rebuilds the `v0.1.0` Figure 1-9 layouts, independent
   Figure S8-S10 outputs, and the eight-page Figure S1-S8 regression package
   from repository-distributed derived tables and approved publication
   intermediates. The submission-asset verifier separately checks the
   ten-page Figure S1-S10 PDF submitted as Additional file 2. Independent
   Figure S9 and S10 files remain regression references, not separate
   Additional files.
4. **Scientific recomputation** reruns analytical models from public,
   author-staged, or licensed inputs. Some workflows require separately
   obtained GEO files, reviewed design metadata, official model assets,
   prepared structural inputs, or licensed database exports.

To verify the current submission package, run:

```powershell
$submissionPackage = $env:PUBLICATION_SUBMISSION_PACKAGE
if (-not $submissionPackage) { throw "Set PUBLICATION_SUBMISSION_PACKAGE to the flat submission-package directory." }
py -3 -B scripts/publication_tables/verify_current_submission_assets.py `
  "$submissionPackage" `
  --inspect-workbooks
```

To rerender the current Figure 8 from the packaged derived tables, run the
following command from the repository root with a new output directory:

```powershell
$figure8Output = Join-Path $env:TEMP ("AA_current_Figure8_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
& (Join-Path $env:ProgramFiles "R\R-4.5.2\bin\Rscript.exe") `
  scripts/molecular_dynamics/render_current_figure8_md.R `
  "$figure8Output"
```

The final 16-page supplementary PDF is held with the submission package. To
apply the reviewed visual-only correction to a separately staged source PDF,
run the R script below with a new output directory. The correction moves the
Figure S6 legend outside the plotting region and removes the redundant page
title from Figure S7; pages 1-5 and 8-16 are checked for pixel equality.

```powershell
$supplementOutput = Join-Path $env:TEMP ("AA_current_supplement_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
& (Join-Path $env:ProgramFiles "R\R-4.5.2\bin\Rscript.exe") `
  scripts/figure_packaging/build_current_supplementary_visual_corrections.R `
  "$env:PUBLICATION_SUPPLEMENTARY_PDF" `
  "$supplementOutput"
```

The archived publication replay remains available on Windows with Python
3.10+ and R 4.5.x:

Install the dedicated replay stack from
`environment/publication_replay_python_20260806.txt`; do not substitute the
Geneformer GPU environment for this rendering-only command.

```powershell
$replayRoot = Join-Path $env:TEMP ("AA_replay_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
$rscript = Join-Path $env:ProgramFiles "R\R-4.5.2\bin\Rscript.exe"
py -3 -B scripts/publication_figures/reproduce_all_publication_figures.py `
  --output-root "$replayRoot" `
  --rscript "$rscript"
```

The output root must not already exist. Figures 1, 2, 4, 5, 6, 8, 9, S9 and
S10 require exact SHA-256 equality. Figures 3, 7 and S8 use documented pixel
tolerances because supported Matplotlib/R rasterizer versions and reference
image encoding can change antialiasing or PNG metadata without changing
source values or panel geometry.

To verify the separately held archived `v0.1.0` publication assets, use the
archive-specific verifier. It accepts either a flat package root containing
the 15 archived assets or the earlier nested layout.

```powershell
$submissionPackage = $env:PUBLICATION_SUBMISSION_PACKAGE
if (-not $submissionPackage) { throw "Set PUBLICATION_SUBMISSION_PACKAGE to the submission-package directory." }
py -3 -B scripts/publication_tables/verify_submission_assets.py `
  "$submissionPackage" `
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
runtime and is not part of the documented runtime validation.

## Geneformer assets and inputs

1. Obtain the public GSE247531 input from NCBI GEO and an author-reviewed design
   table with the columns documented in `config/input_paths.example.tsv`.
2. Run the balanced-input preparation script with a new output directory:

```bash
Rscript scripts/geneformer/70_aa_geneformer_prepare_balanced_cd34_public_20260804_v1.R \
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

The current manuscript uses Figure 1-8. Exact current references are stored in
`reference_outputs/current_manuscript/`; source and access boundaries are in
`CURRENT_MANUSCRIPT_FIGURE_SOURCE_MAP.tsv`. Figure 8 is rerenderable from the
current five-complex derived tables and reproduces the submitted geometry to
the documented pixel tolerance. Figure 7 is retained as an exact final-image
reference because prepared receptor, ligand, pose and visualization assets are
not redistributed. Supplementary Figures S1-S16 are verified as a separately
held 16-page PDF using the current submission-asset verifier.

The archived `v0.1.0` references remain available for historical regression:

The archived `v0.1.0` publication appearance of Figure 1-9 and Figure S1-S10
is replayable or regression-verifiable. This does **not** mean that every analysis
is self-contained from raw data. Figures 1, 2, 4, 5 and 6 and several legacy
supplementary panels still require author-staged or licensed upstream inputs
for scientific recomputation. Archived Figure 7 is rerendered from the complete derived
summary of five 100 ns systems (500 ns total), while raw trajectories remain
excluded. Figure 9/S9/S10 is rerendered from recorded model outputs; repeating
model inference requires official Geneformer assets. Exact statuses are in
`REPRODUCIBILITY_MATRIX.tsv` and `FIGURE_SOURCE_MAP.tsv`.

The archived Figure 9/S9/S10 presentation pipeline separates panel rendering from the
final RGB conversion. This preserves black titles and a white background
without recomputing analytical results.

## License and citation

Repository-authored code is BSD-3-Clause. Third-party datasets, databases,
Geneformer code/model assets, SPSS, Cytoscape, GROMACS, AutoDock Vina and
R/Python packages retain their own terms. The repository URL and issued Zenodo
DOI are recorded in `CITATION.cff`.

Software availability information for manuscript reporting:

- Project name: `aa-donor-recovery-perturbation`.
- Project home page: <https://github.com/yuyi-wang-imu/aa-donor-recovery-perturbation>.
- Archived version: `v0.1.0`, <https://doi.org/10.5281/zenodo.21837457>.
- Operating systems: Windows for the byte-stable publication replay, and
  Windows/Linux for analytical scripts where dependencies permit.
- Programming languages: Python 3.10+ and R 4.5.x.
- Other requirements: workflow-specific dependencies and external inputs are
  documented in `environment/`, `config/input_paths.example.tsv`, and
  `DATA_AND_LICENSES.md`.
- License: BSD-3-Clause for repository-authored code.
- Restrictions: external software, data, models, licensed databases, prepared
  structures and raw molecular-dynamics trajectories retain their own access
  terms and are not redistributed here.

The DOI identifies the archived `v0.1.0` release. Subsequent maintenance on
the default branch does not alter that archived release or its DOI. The current
manuscript title and evidence hierarchy above are a documentation alignment;
they do not silently replace the contents of the archived release.
