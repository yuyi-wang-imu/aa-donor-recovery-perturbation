# BMC Pharmacology and Toxicology simulated-editor v5 alignment

This maintenance candidate aligns public repository documentation with the
simulated-editor v5 manuscript without moving the immutable v0.1.0 tag or
changing DOI 10.5281/zenodo.21837457.

## Scientific wording aligned

- Replaces `prescription-derived compounds` with database compounds associated
  with herbs recurrent in published aplastic-anemia prescriptions.
- States that those compounds were not chemically quantified or verified in an
  administered formulation.
- Distinguishes 18 compounds with retained SwissTargetPrediction records from
  23 unique structurally tractable docking ligands.
- Describes the donor coordinate as a follow-up-associated healthy-directed
  longitudinal transcriptomic trajectory rather than a validated clinical
  recovery endpoint.
- Records that only TOP2A and GSK3B had donor-bootstrap intervals excluding
  zero, while neither passed matched-background multiplicity correction and the
  pooled matched-null test was nonsignificant.
- Keeps docking, molecular dynamics and scTenifoldKnk as downstream,
  hypothesis-generating analyses rather than candidate-selection evidence.

## Archive boundary

- Zenodo v0.1.0 remains the immutable code and selected-derived-data archive at
  commit `b850d7e6a6141962d3a2d38021d0780e3c5907a1`.
- The old DOI does not contain this maintenance candidate or the complete
  13.2 GB MD production-record inventory.
- Current repository documentation and final package checksums will be tied to
  the new commit SHA only after the v5 package passes final visual and
  structural QA.

## Zenodo metadata-only description candidate

The public record may use the following description without changing the
archived v0.1.0 ZIP, tag, DOI or version:

> Reproducible analysis code and selected derived source tables supporting a donor-aware longitudinal transcriptomic study of aplastic anemia. The workflow integrates traceable candidate provenance from structurally specified database compounds associated with herbs recurrent in published aplastic-anemia prescriptions, CD34+ hematopoietic and bone-marrow single-cell context, a healthy-directed longitudinal transcriptomic trajectory, and paired in silico deletion and overexpression analyses using Geneformer. The database compounds were not chemically quantified or verified in an administered formulation. scTenifoldKnk, molecular docking and molecular-dynamics simulations are complementary downstream, hypothesis-generating analyses. Release v0.1.0 corresponds to Git commit b850d7e6a6141962d3a2d38021d0780e3c5907a1 and is distributed under the BSD-3-Clause license. Later GitHub states, including the BMC Pharmacology and Toxicology simulated-editor v5 manuscript alignment, are not contained in the Zenodo v0.1.0 archive. Public transcriptomic inputs are available from NCBI GEO under GSE247531, GSE165870 and GSE145668. Database-licensed records, third-party software and the 13.2 GB complete molecular-dynamics production-record inventory are not redistributed in this archive; derived time series, summaries and analysis scripts are public in the linked GitHub repository.

Retain the existing title, creator `Wang, Yuyi`, ORCID
`0009-0009-4087-7142`, resource type `Software`, version `v0.1.0`,
BSD-3-Clause licence, GitHub relation and DOI. After saving metadata, verify
that the archived file remains
`yuyi-wang-imu/aa-donor-recovery-perturbation-v0.1.0.zip`, size 48,051,582
bytes and MD5 `c54383bbd520f980095dedfe5334c156`.
