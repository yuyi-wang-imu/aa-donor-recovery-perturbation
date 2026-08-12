# Data and software redistribution boundaries

Repository-authored code is released under BSD-3-Clause. This license does not
relicense external data, database exports, software, pretrained models, or
author-controlled clinical/design metadata.

## Included

- Repository-authored scripts, configuration templates, and validation tools.
- Selected derived numerical tables required to inspect or rerender manuscript
  figures, where the table does not reproduce a restricted raw database export.
- Aggregated and model-derived outputs for the Geneformer calibration panels.
- Publication-facing image/PDF intermediates and reference figures used only
  for layout replay and regression testing. These assets are
  not treated as substitutes for analytical source data.
- Derived public-study tables may retain pseudonymous subject labels (for
  example, UPN codes) and single-cell barcodes solely for traceable provenance;
  they contain no names, contact details, or private design table.

## Obtain separately

- NCBI GEO data: GSE247531, GSE165870, and GSE145668. Users must retrieve public
  files from GEO and follow NCBI terms and the originating study's conditions.
- Geneformer source and Geneformer-V2-104M assets from the official upstream
  distribution. The upstream model card labels the model Apache-2.0. The
  approximately 418 MB `model.safetensors` file is excluded from this repository.
- R/Python packages, Cytoscape, GROMACS, AutoDock Vina, and SPSS under their own
  licenses. This repository does not distribute those programs.

Recorded Geneformer provenance:

- model: `Geneformer-V2-104M`
- model SHA-256: `fff5cba29ddd8792991fa77b4872246fbe548a178cebda3775cdc72b67780e7f`
- official-source snapshot commit recorded by the analysis: `04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5`
- upstream model page: `https://huggingface.co/ctheodoris/Geneformer`

The public repository must link to upstream assets rather than commit model
weights, dictionary snapshots, or a vendored Geneformer checkout.

## Not redistributed

- The complete publication supplementary workbooks are not duplicated into GitHub.
  `PUBLICATION_ASSET_CHECKSUMS.tsv` and
  `scripts/publication_tables/verify_submission_assets.py` verify the separate
  journal files, including 390 prescription records and the corrected 12-rule
  table, without republishing their archived database rows as source exports.
- TCMSP, SwissTargetPrediction, OMIM, GeneCards, TTD, DisGeNET, bibliographic,
  or other licensed/raw database exports. Derived associations must be reviewed
  case by case and cannot be treated as permission to republish the source dump.
- Private participant-level design metadata. Users must reconstruct required
  metadata from the cited public GEO studies and review it under applicable
  data-use and privacy requirements.
- Prepared receptor/ligand structures and exact pose-visualization project files.
- Raw molecular-dynamics topology, checkpoint, energy and trajectory files
  (`*.tpr`, `*.cpt`, `*.edr`, `*.xtc`, `*.trr`, `*.dcd`, `*.nc`).
- Secrets, credentials, API tokens, private keys, local `.env` files, personal
  paths, model weights, caches, and generated output directories.

## Reproducibility consequence

The repository provides a complete publication-replay layer for the submitted
figure set and a scientific-recomputation layer wherever legal, public, or
author-reviewed inputs can be staged. It is not a self-contained from-raw
capsule for every analysis. Author-prepared structures, licensed source
exports, raw MD trajectories, official model assets, and reviewed participant
design staging remain external. The exact boundary for every figure and table
is recorded in `REPRODUCIBILITY_MATRIX.tsv` and must remain visible in the
public README.
