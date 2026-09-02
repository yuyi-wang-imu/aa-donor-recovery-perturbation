# BMC simulated-editor v5 data and code availability

## Ready-to-paste manuscript wording

### Availability of data and materials

The curated prescription dataset and processed outputs are provided in
Additional files 2, 3 and 5-6; the search and eligibility framework is provided
in Additional file 7. Public transcriptomic data were obtained from the NCBI
Gene Expression Omnibus under GSE247531, GSE165870 and GSE145668; no new
primary human transcriptomic data were generated. Analysis code and selected
derived source tables, including data underlying Fig. 8 and Supplementary
Figs. S14-S16, are available at
https://github.com/yuyi-wang-imu/aa-donor-recovery-perturbation. The immutable
software and derived-data release v0.1.0 remains archived at Zenodo
(https://doi.org/10.5281/zenodo.21837457; commit
`b850d7e6a6141962d3a2d38021d0780e3c5907a1`); it does not by itself archive
subsequent manuscript-specific repository revisions, which are identified by
commit SHA in the repository release notes. Repository-authored code is
distributed under the BSD-3-Clause license. External software, model assets and
licensed databases retain their original access terms. Complete saved
production records for five ligand-containing 100 ns trajectories and one
matched 100 ns HIF1A-ARNT apo reference (six trajectories; 600 ns total;
13.2 GB) are not included in the Additional files or the archived v0.1.0 ZIP
because of their aggregate size. Editors and peer reviewers may request the
complete retained files through the institutional corresponding-author route
(Wendurige, wendurige@imu.edu.cn) for verification of the reported analyses;
requests should identify the manuscript and requested systems, and transfer
will use an institutionally approved secure method appropriate to the file
size. Derived time series, final-window summaries, analysis scripts and
available simulation-input metadata are public in the repository.

## Repository and citation actions before public synchronization

- Freeze the final v5 manuscript, cover letter and any modified Additional file.
- Refresh the 17-file BMC package checksum table from those exact files.
- Rebuild `MANIFEST.tsv` after every repository document and checksum update.
- Run the repository structural validators and record the exact commit SHA.
- Cite Zenodo v0.1.0 only as the immutable software and selected-derived-data
  release; cite the current GitHub commit separately for v5 documentation.
- Do not create or cite a new MD DOI unless the complete intended record has
  actually been deposited and published.

## Current access boundary

- Public now: public GEO accession routes, repository-authored code, selected
  derived source tables, Figure 8 and Supplementary Figs. S14-S16 replay data,
  summaries, manifests and available input metadata.
- Not public in Zenodo v0.1.0: later manuscript-specific commits and the 13.2 GB
  complete six-system MD production records.
- Request route: editors and peer reviewers contact the institutional
  corresponding author and identify the manuscript plus requested systems;
  transfer uses an institutionally approved secure method appropriate to the
  data volume.
