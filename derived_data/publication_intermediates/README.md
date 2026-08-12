# Publication-facing intermediate assets

This directory contains author-generated intermediate PNGs and small derived
tables used to replay the final manuscript figure layouts. These files are
not raw clinical data, licensed database exports, prepared docking structures,
model weights, or molecular-dynamics trajectories.

The exact-final replay layer is deliberately separate from the analytical
layer. Analytical scripts regenerate scientific tables from public or
author-staged inputs. The replay layer then applies only approved composition,
spacing, color-mode, and panel-letter operations so the published figures can
be checked against `reference_outputs/`.

Files named `Figure*_pre_*` are approved publication intermediates. The CSV
files contain only figure-level derived node, edge, or annotation values that
are already represented in the manuscript figures or supplementary tables.
They do not contain the full TCMSP or SwissTargetPrediction source exports.

`supplementary_pages/` contains the eight accepted, single-page PDF assets
used by the packaging-only Supplementary Figures S1-S8 workflow. They support
publication-package replay; they are not represented as raw analytical inputs.
