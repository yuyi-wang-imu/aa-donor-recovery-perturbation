# GitHub and Zenodo release checklist

This checklist is for a new repository. It does not authorize deletion,
unpublishing, overwriting or DOI withdrawal for the historical repository.

## Repository preparation

- [x] Repository created at
      `https://github.com/yuyi-wang-imu/aa-donor-recovery-perturbation`.
- [x] Figure 9/S9/S10 use black titles, white backgrounds and RGB output.
- [x] Private participant design metadata excluded; included public-study
      pseudonymous labels documented in `DATA_AND_LICENSES.md`.
- [x] Repository title, URL and currently confirmed software author recorded in
      `CITATION.cff` without inventing a DOI.
- [x] Add publication intermediates and deterministic replay for
      Figures 1, 4, 5 and 6 while retaining conditional upstream-analysis labels.
- [x] Confirm that TCMSP, SwissTargetPrediction and other licensed/raw database
      exports are absent and that only permissible derived tables remain.
- [x] Rebuild Figure 1-9 and Figure S8-S10 and verify exact checksums or documented
      pixel tolerances against the publication submission references.
- [x] Confirm the submitted Figure 5 and Figure 8 reference images and record
      the Additional file 4 checksum without redistributing the journal workbook.

## Local static checks

```bash
python3 scripts/build_manifest.py
python3 scripts/run_workflow.py --check
python3 scripts/validate_repository.py
```

The checks are non-analytical. They do not rerun MD, Geneformer, WGCNA,
scTenifoldKnk or docking. R parse validation must additionally be run on a
machine with `Rscript` available. On Windows, use the installed
absolute executable path if `Rscript` is not on `PATH`.

## Repository publication and archiving

- [x] Add only paths listed in `MANIFEST.tsv`.
- [x] Keep the repository private for a final GitHub secret/size/license review.
- [x] Add the confirmed new repository URL to `CITATION.cff` and README.
- [x] Run a clean-clone structural check and the non-analytical Python/R checks.
- [ ] Where external inputs are available, run the
      minimal smoke tests in a new output directory.
- [x] Push and verify the final approved commit while the repository remains private.
- [x] Author approval received for public visibility after repository checks pass.
- [x] Enable the new repository in Zenodo before creating its first formal release.
- [x] Create tag `v0.1.0`, matching `CITATION.cff`, on the verified final commit;
      do not reuse any historical tag.
- [x] Create the GitHub Release and wait for Zenodo to archive that release and mint a new DOI.
- [x] Add the issued DOI only after it exists; never guess a DOI.
- [x] Rebuild `MANIFEST.tsv` and the release archive after metadata changes.

## Explicit prohibitions

- Do not force-push. Change repository visibility only after author approval and
  successful repository checks.
- Do not delete, archive, or change visibility of the old repository without a
  separate, explicit author decision.
- Do not withdraw, replace or repoint the historical DOI.
- Do not upload raw MD trajectories, model weights, secrets, private metadata,
  licensed source exports, local caches, or the surrounding manuscript folders.
