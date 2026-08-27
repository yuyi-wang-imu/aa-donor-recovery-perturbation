# Companion molecular-dynamics dataset

The complete saved production records underlying Figure 8 and Supplementary
Figures S14-S16 are archived as an independent Zenodo dataset:

<https://doi.org/10.5281/zenodo.22131869>

The existing software DOI <https://doi.org/10.5281/zenodo.21837457> remains the
immutable `v0.1.0` code and derived-data snapshot and is not replaced by the MD
dataset DOI.

## Dataset scope

The companion record contains six ZIP64 system packages:

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

An additional archive contains the exact final analysis-ready source tables for
Figure 8 and Supplementary Figures S14-S16. The GitHub copies under
`derived_data/molecular_dynamics/` remain the lightweight figure-replay layer;
the Zenodo record is the immutable raw-system layer.

These simulations describe conformational behavior under the specified
computational conditions. They do not establish biochemical affinity, target
engagement, inhibition or therapeutic efficacy. The HIF1A-ARNT apo/butin
comparison is exploratory and uses the engineered PDB 4H6J PAS-B heterodimer.
