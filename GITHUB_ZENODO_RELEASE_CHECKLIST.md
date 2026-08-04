# New GitHub/Zenodo release checklist

This checklist is for a new repository. It does not authorize deletion,
unpublishing, overwriting or DOI withdrawal for the historical repository.

## Before the author creates the repository

- [ ] Resolve the Figure 9/S9/S10 v5 versus later v6/v7 presentation-version conflict.
- [ ] Review participant/design metadata before any metadata file is added.
- [ ] Confirm author list, software title and intended new repository name.
- [ ] Decide whether Figures 1, 4, 5 and 6 will remain documented as conditional
      reproducibility or whether their missing accepted inputs/renderers will be added.
- [ ] Confirm that TCMSP, SwissTargetPrediction and other licensed/raw database
      exports are absent and that only permissible derived tables remain.

## Local static checks

```bash
python scripts/build_manifest.py
python scripts/run_workflow.py --check
python scripts/validate_repository.py
```

The checks are non-analytical. They do not rerun MD, Geneformer, WGCNA,
scTenifoldKnk or docking. R parse validation must additionally be run on a
machine with `Rscript` available; it is not available in the current
local environment.

## After the author creates the new private repository

- [ ] Add only paths approved by `PROPOSED_UPLOAD_MANIFEST_20260804_v1.tsv`.
- [ ] Keep the repository private for a final GitHub secret/size/license review.
- [ ] Add the confirmed new repository URL to `CITATION.cff` and README.
- [ ] Run a clean-clone structural check and, where inputs are available, the
      minimal smoke tests in a new output directory.
- [ ] Make public only after manuscript/code-version approval.
- [ ] Create a new release tag chosen by the author; do not reuse the old release tag.
- [ ] Enable Zenodo only for the new repository and mint a new DOI.
- [ ] Add the issued DOI only after it exists; never guess a DOI.
- [ ] Rebuild `MANIFEST.tsv` and the release archive after metadata changes.

## Explicit prohibitions

- Do not push or publish from this audit thread.
- Do not delete or hide the old repository.
- Do not withdraw, replace or repoint the historical DOI.
- Do not upload raw MD trajectories, model weights, secrets, private metadata,
  licensed source exports, local caches, or the surrounding manuscript folders.
