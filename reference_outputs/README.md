# Golden publication outputs

These files are byte-level or pixel-level references extracted from the frozen
manuscript submission package. They are used only for regression testing of
the public figure-reproduction workflow. They are not treated as analytical
inputs and must never substitute for the source tables and scripts.

Most regenerated figures are expected to match byte-for-byte. Figures 3, 7,
and S8 are checked by documented pixel tolerances because frozen re-encoding
or supported Matplotlib/R rasterizer versions change antialiasing or PNG
metadata without changing the source values or panel geometry.
