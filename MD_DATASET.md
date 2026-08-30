# Molecular-dynamics production-record inventory and access

The complete saved production records underlying Figure 8 and Supplementary
Figures S14-S16 are not stored in GitHub because the retained six-system
inventory is 13.2 GB. They are available on reasonable request from Wendurige
(`wendurige@imu.edu.cn`) and will be provided to editors and reviewers during
peer review. No separate public DOI has been issued for these records.

The existing software DOI <https://doi.org/10.5281/zenodo.21837457> remains the
immutable `v0.1.0` code and derived-data snapshot.

## Dataset scope

The retained inventory contains six system packages:

1. TOP2A-sesamin, 100 ns;
2. GSK3B-linarin, original 0-50 ns and 50-100 ns segments;
3. KIT-3'-O-methylorobol, 100 ns;
4. HIF1A-ARNT-butin, exploratory 100 ns trajectory;
5. SYK-isofucosterol, original 0-50 ns and 50-100 ns segments; and
6. a matched HIF1A-ARNT apo reference, 100 ns.

Together these records represent 600 ns of production simulation. Every system
package contains the complete saved trajectory at the recorded 10 ps interval,
the matching production TPR, energy, checkpoint and final-coordinate files,
run-specific topologies and parameters, prepared starting structures, protocol
files, path-sanitized production logs, file-level provenance, licences,
manifests and SHA-256 checksums. The PubChem source SDF downloads are not
redistributed; their CIDs, source URLs and source-file SHA-256 values are
recorded instead.

The exact final analysis-ready source tables for Figure 8 and Supplementary
Figures S14-S16 are distributed under `derived_data/molecular_dynamics/` as the
lightweight figure-replay layer. The request-based records are the complete
saved-system layer.

These simulations describe conformational behavior under the specified
computational conditions. They do not establish biochemical affinity, target
engagement, inhibition or therapeutic efficacy. The HIF1A-ARNT apo/butin
comparison is exploratory and uses the engineered PDB 4H6J PAS-B heterodimer.
