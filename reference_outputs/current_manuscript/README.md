# Current-manuscript image references

This directory contains the eight visually approved, lowercase-panel main-
figure regression references from the final Scientific Reports transfer
package. The submission-file byte sizes and SHA-256 values and the separately
frozen regression-reference byte sizes and SHA-256 values are both recorded in
`CURRENT_MANUSCRIPT_ASSET_CHECKSUMS.tsv`. The v4 submission PNGs are
pixel-identical to these references but carry standardized 300 dpi metadata;
the frozen copies retain their previously validated PNG encoding.

The photographic graphical abstract is a retained asset from the prior Human
Genomics submission and is **not** part of the Scientific Reports transfer
package. Its image sources and reuse boundaries are documented in
`GRAPHICAL_ABSTRACT_IMAGE_PROVENANCE.md`.

These images are regression references, not analytical inputs. Figure 8 can be
rerendered from the repository-distributed derived tables with
`scripts/molecular_dynamics/render_current_figure8_md.R`. Figure 7 is retained
as a checksum-verifiable final reference because the prepared receptor,
ligand, pose and molecular-visualization assets are intentionally not
redistributed. The source and access boundaries for every current figure are
specified in `CURRENT_MANUSCRIPT_FIGURE_SOURCE_MAP.tsv` and
`CURRENT_MANUSCRIPT_REPRODUCIBILITY_MATRIX.tsv`.

The final v12 manuscript contains 104 sequentially cited references;
manuscript-only revisions did not change the eight approved figure binaries or
the archived analysis release. The archived `v0.1.0` Figure 1-9
references remain unchanged in the sibling
`main_figures` directory. Do not use the archived numbering to verify the
current eight-figure manuscript.
