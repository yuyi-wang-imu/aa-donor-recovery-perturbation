# Data and software redistribution boundaries

Repository-authored code is released under BSD-3-Clause. This license does not
relicense external data, database exports, software, pretrained models, or
author-controlled clinical/design metadata.

## Included

- Repository-authored scripts, configuration templates, and validation tools.
- Selected derived numerical tables required to inspect or rerender manuscript
  figures, where the table does not reproduce a restricted raw database export.
- Aggregated and model-derived outputs for the Geneformer calibration panels.

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

- Author-curated 390-record prescription workbook and verified SPSS exports.
- TCMSP, SwissTargetPrediction, OMIM, GeneCards, TTD, DisGeNET, bibliographic,
  or other licensed/raw database exports. Derived associations must be reviewed
  case by case and cannot be treated as permission to republish the source dump.
- Participant-level design metadata until the author confirms that only public
  GEO identifiers and non-sensitive fields remain.
- Prepared receptor/ligand structures and exact pose-visualization project files.
- Raw molecular-dynamics topology, checkpoint, energy and trajectory files
  (`*.tpr`, `*.cpt`, `*.edr`, `*.xtc`, `*.trr`, `*.dcd`, `*.nc`).
- Secrets, credentials, API tokens, private keys, local `.env` files, personal
  paths, model weights, caches, and generated output directories.

## Reproducibility consequence

The repository supports code and selected-derived-data verification. It is not
yet a complete from-raw public capsule for every panel: Figures 1, 4, 5 and 6
and parts of the legacy supplementary package require author-controlled or
third-party inputs. These limitations are listed in `FIGURE_SOURCE_MAP.tsv` and
must remain visible in the public README.
