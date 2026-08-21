# Current-manuscript image references

This directory contains the eight main figures and the graphical abstract used
by the current Human Genomics submission package. Their filenames, byte sizes
and SHA-256 values are recorded in `CURRENT_MANUSCRIPT_ASSET_CHECKSUMS.tsv`.

These images are regression references, not analytical inputs. Figure 8 can be
rerendered from the repository-distributed derived tables with
`scripts/molecular_dynamics/render_current_figure8_md.R`. Figure 7 is retained
as a checksum-verifiable final reference because the prepared receptor,
ligand, pose and molecular-visualization assets are intentionally not
redistributed. The source and access boundaries for every current figure are
specified in `CURRENT_MANUSCRIPT_FIGURE_SOURCE_MAP.tsv` and
`CURRENT_MANUSCRIPT_REPRODUCIBILITY_MATRIX.tsv`.

The archived `v0.1.0` Figure 1-9 references remain unchanged in the sibling
`main_figures` directory. Do not use the archived numbering to verify the
current eight-figure manuscript.
