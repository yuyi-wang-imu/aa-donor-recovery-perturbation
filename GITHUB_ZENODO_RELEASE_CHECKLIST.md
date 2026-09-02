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

## Post-release metadata maintenance (2026-08-19)

- [x] Keep tag `v0.1.0` fixed at commit
      `b850d7e6a6141962d3a2d38021d0780e3c5907a1`.
- [x] Keep version DOI `10.5281/zenodo.21837457`; do not create a replacement
      release or silently repoint the DOI.
- [x] Align default-branch README and citation metadata with the current
      Human Genomics manuscript title and evidence hierarchy.
- [x] State explicitly that the archived `v0.1.0` replay retains its legacy
      Figure 1-9/Figure S1-S10 numbering, whereas the current manuscript uses
      eight main figures and Figure S1-S16.
- [x] If Zenodo metadata is edited, change only title, description, and
      keywords, then verify the public record still resolves to the same DOI
      and the same archived `v0.1.0` file.
- [x] Do not replace archived files or create a new version merely to make a
      documentation-only maintenance commit appear inside the old ZIP.

Verified on 2026-08-20: the public record retained version DOI
`10.5281/zenodo.21837457`, concept DOI `10.5281/zenodo.21837456`, version
`v0.1.0`, and archived file
`yuyi-wang-imu/aa-donor-recovery-perturbation-v0.1.0.zip` (48,051,582 bytes;
MD5 `c54383bbd520f980095dedfe5334c156`).

## Current Human Genomics manuscript alignment (2026-08-21)

- [x] Preserve the archived Figure 1-9 and Figure S1-S10 replay without
      moving tag `v0.1.0` or changing DOI `10.5281/zenodo.21837457`.
- [x] Add exact current references for Figure 1-8 and the graphical abstract.
- [x] Add a current 19-file submission checksum table and a dedicated verifier.
- [x] Add current Figure 1-8 source and reproducibility maps and the Figure
      S1-S16 page map.
- [x] Add the five-candidate Figure 8 derived source tables and a path-neutral
      R renderer; do not add raw trajectories, topologies or checkpoints.
- [x] Verify the separately held current submission package without copying the
      manuscript, cover letter or supplementary workbooks into GitHub.
- [x] Keep prepared docking structures, ligand/pose projects, licensed source
      exports, model weights and private design metadata outside the repository.

## Scientific Reports transfer alignment (2026-08-26)

- [x] Create the local branch `scientific-reports-transfer-20260826` from
      `74d91070c762bd464c104d57d38ee1eb97a1e70a` without changing remote state.
- [x] Preserve the immutable `v0.1.0` tag, archived ZIP and version DOI
      `10.5281/zenodo.21837457`; do not create a new tag, release or DOI for this
      journal transfer.
- [x] Retain the Human Genomics package history as provenance while defining a
      separate 16-file Scientific Reports transfer package.
- [x] Exclude the prior graphical abstract and standalone main-table files from
      the Scientific Reports package; retain their scientific content or
      provenance records where applicable.
- [x] Define strict first-citation renumbering for Supplementary Figs. S1-S16
      and Supplementary Tables S1-S10 in the current mapping files.
- [x] Correct the retained Cuscutae Semen image provenance to Song Yang et al.
      (2016), DOI `10.1155/2016/8656740`, Figure 1 crude product.
- [x] Replace the eight current main-figure regression references with the
      visually approved lowercase-panel Scientific Reports files from the
      final package. The v4 submission PNGs are pixel-identical to these frozen
      references and carry the final 300 dpi metadata; the frozen repository
      copies retain their validated original PNG encoding.
- [x] Replace `CURRENT_MANUSCRIPT_ASSET_CHECKSUMS.tsv` with the exact final v12
      16-file Scientific Reports inventory after all files passed review.
- [x] Run `verify_scientific_reports_submission_assets.py --inspect-workbooks`
      against a clean copy of the final package and confirm strict
      supplementary first-use order, absence of `Additional file` terminology
      and absence of forbidden graphical-abstract/main-table extras.
- [x] Rebuild `MANIFEST.tsv` last, run all repository validators, review the
      staged diff, and confirm author authorization before commit and push.

The final v12 checksum mapping was prepared and validated locally on 2026-08-28.
The repository's frozen publication-reference binaries were intentionally left
unchanged because the v4 figures are pixel-identical and differ only in PNG
metadata. The author explicitly authorized the Scientific Reports repository
update on 2026-08-26. The archived Zenodo version DOI, files and tag remain
unchanged; only the existing record metadata may identify the live transfer
branch and final manuscript-asset mapping.

## Simulated-editor v5 BMC alignment (2026-09-02)

- [x] Preserve the immutable Zenodo v0.1.0 tag, ZIP, version DOI
      `10.5281/zenodo.21837457`, and archived commit `b850d7e6...`.
- [x] State that the old DOI does not contain later manuscript-specific commits
      or the 13.2 GB molecular-dynamics production-record inventory.
- [x] Clarify the 18 SwissTargetPrediction compounds versus 23 structurally
      tractable docking ligands and state that the database compounds were not
      chemically verified in an administered formulation.
- [x] Record the matched-background negative calibration results and retain the
      five genes as exploratory directional candidates.
- [x] Document the editor/reviewer request route for complete MD production
      records without claiming that those files are public or archived.
- [ ] Freeze the final v5 manuscript, cover letter and any revised Additional
      file after visual QA.
- [ ] Regenerate `BMC_PHARMACOLOGY_TOXICOLOGY_ASSET_CHECKSUMS.tsv` from the
      frozen v5 package and rebuild `MANIFEST.tsv` last.
- [ ] Run `scripts/run_workflow.py --check` and
      `scripts/validate_repository.py` on the staged v5 branch.
- [ ] Review the exact Git diff and obtain author confirmation of the commit
      hash before pushing or merging to the public default branch.
- [ ] If the Zenodo public metadata is updated, change metadata only and verify
      that DOI `10.5281/zenodo.21837457`, version v0.1.0 and its archived ZIP
      remain unchanged. Do not claim that metadata editing changes archive
      contents.
