# Current Figure 8 derived source data

These three files support the current BMC Pharmacology and Toxicology Figure 8 without
placing multi-gigabyte molecular-dynamics trajectories, topologies, checkpoints
or server logs in GitHub. Six complete saved 0-100 ns system records (13.2 GB;
600 ns total) are available on reasonable request from Wendurige
(`wendurige@imu.edu.cn`) and to editors and reviewers during peer review. No
separate public DOI has been issued for these records.

- `Figure8_five_candidates_time_series_source.tsv.gz` contains 10,001 frames
  from 0 to 100 ns for each of five ligand-containing complexes and four
  metrics. `time_ns` is simulation time, `value` is the frame-level metric and
  `smooth` is the centered 1.01 ns rolling median used for display.
- `Figure8_five_candidates_ca_rmsf_source.tsv` contains C-alpha RMSF values by
  residue index. The chain column identifies the HIF1A and ARNT segments in the
  heterodimer; residue indices are system-specific and are not homologous
  coordinates across proteins.
- `Figure8_final20ns_quantitative_summary.tsv` contains the reported 80-100 ns
  summaries: mean protein-backbone RMSD, mean ligand heavy-atom RMSD after
  protein fitting, median minimum protein-ligand heavy-atom distance, mean
  ligand proximity fraction and median C-alpha RMSF.

The five complexes are TOP2A-sesamin, GSK3B-linarin,
KIT-3'-O-methylorobol, HIF1A-ARNT-butin and SYK-isofucosterol. The
HIF1A-ARNT-butin trajectory is an exploratory comparison from a defined
starting structure. Figure 8 describes conformational behavior during the
sampled trajectories; it is not a comparative binding-affinity estimate.

Run `scripts/molecular_dynamics/render_current_figure8_md.R` from the
repository root with a new output directory. The renderer validates case
membership, frame counts, time range, finite values and proximity bounds before
creating the figure.
