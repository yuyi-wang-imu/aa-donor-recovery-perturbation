# Supplementary Figures S14-S16 derived source data

This directory contains the minimal derived numeric source tables used for
Supplementary Figures S14-S16 in the Scientific Reports manuscript. These are
panel-level analysis outputs, not raw trajectories. Multi-gigabyte
molecular-dynamics trajectories, topologies, checkpoints and production logs
are not committed to GitHub; the complete saved 0-100 ns records are archived
in the companion Zenodo dataset at
<https://doi.org/10.5281/zenodo.22131869>.

Frame-level time series span 0-100 ns at 0.01-ns (10-ps) intervals. A complete
single-system series therefore contains 10,001 frames. The paired S16 time-series
table contains 20,002 rows: 10,001 for the apo system and 10,001 for the
butin-starting-pose system. The late window used for the reported final-20-ns
summaries and RMSF profiles is 80-100 ns, inclusive.

## Supplementary Figure S14: TOP2A quality control

Files are under `S14_TOP2A_QC/`.

- Panel a: `coordinate_metrics_timeseries.tsv`, column
  `protein_radius_of_gyration_nm`.
- Panel b: `sasa.tsv`, column `protein_sasa_nm2`.
- Panel c: `temperature_pressure_density.tsv`, column `temperature_K`.
- Panel d: `temperature_pressure_density.tsv`, column `density_kg_m3`.
- Panel e: `temperature_pressure_density.tsv`, column `pressure_bar`.
- Panel f: `residue_contact_occupancy.tsv`. The ten displayed pocket residues
  were selected by descending maximum occupancy across the full 0-100-ns and
  final 80-100-ns windows, with `new_resid` used to break ties. Residue labels
  use `original_chain`, `original_resname`, and `original_resid`.

## Supplementary Figure S15: paired HIF1A-ARNT trajectories

Files are under `S15_HIF1A_ARNT_paired/`.

- Panels a-f use the apo and butin-starting-pose
  `all_metrics_timeseries.tsv` tables. The plotted columns are
  `whole_protein_backbone_rmsd_nm`,
  `chain_A_backbone_rmsd_global_fit_nm`,
  `chain_B_backbone_rmsd_global_fit_nm`,
  `interface_residue_contact_pairs_le_0_40_nm`,
  `chain_A_B_com_distance_nm`,
  `interface_buried_sasa_per_side_nm2`, and
  `interface_hydrogen_bond_count`.
- Panel g uses the apo and butin-starting-pose `ca_rmsf_by_chain.tsv` tables,
  column `last_80_100ns_ca_rmsf_global_fit_nm`.
- Panel h uses the butin-starting-pose `all_metrics_timeseries.tsv` table,
  columns `butin_heavy_rmsd_after_protein_fit_nm`,
  `pocket_butin_com_distance_nm`, and `pocket_occupied_le_0_50_nm`.

## Supplementary Figure S16: HIF1A pocket geometry

Files are under `S16_HIF1A_ARNT_pocket/`.

- Panel a: the paired pocket-geometry time-series table, column
  `pocket_heavy_geometric_center_displacement_from_initial_nm`.
- Panel b: the same table, column
  `hif1a_cys255_cys337_sg_distance_nm`.
- Panel c: the butin rows of the same table, columns
  `butin_com_to_initial_pocket_geometric_center_nm` and
  `butin_com_to_dynamic_pocket_geometric_center_nm`.
- Panel d: the pocket-residue RMSF table, column
  `last_80_100ns_ca_rmsf_global_fit_nm`.

The S16 files are LF-normalized, cross-platform serializations of the final
figure-source tables. Their fields, row order, missing values, and numeric
content are unchanged. R serialization represents some values more compactly
(for example, `0` instead of `0.0`) and writes missing ligand quantities as
`NA`; numeric agreement with the analysis tables was verified to within
5.2 x 10^-15. The companion Zenodo analysis-source archive preserves the
original final exports together with file-level checksums.

## Interpretation boundary

Each condition is represented by one molecular-dynamics trajectory. Time frames
are correlated observations, not independent biological or simulation
replicates, and the figures do not support frame-level inferential testing.
These analyses describe conformational behavior conditional on the specified
starting structures; they do not estimate comparative binding affinity or
establish experimental target engagement, inhibition, or therapeutic efficacy.
The HIF1A-ARNT-butin trajectory is an exploratory starting-pose comparison. Its
pocket-center displacement is not a pocket-volume measurement, and the
Cys255-Cys337 sulfur distance is a representative geometric coordinate rather
than a complete measure of pocket opening or closure.
