#!/usr/bin/env Rscript

cran <- c(
  "WGCNA", "dynamicTreeCut", "fastcluster", "Matrix", "ggplot2", "data.table",
  "pheatmap", "ggrepel", "readr", "dplyr", "patchwork", "scales",
  "svglite", "ragg", "tidyr", "jsonlite", "png", "future", "digest",
  "magick", "officer", "Seurat", "scTenifoldNet", "scTenifoldKnk"
)
bioc <- c("limma", "DESeq2")

missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

required <- c(
  cran, bioc
)
still_missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(still_missing)) {
  stop(
    "Required R packages remain unavailable after installation: ",
    paste(still_missing, collapse = ", ")
  )
}

message("PASS: all declared R dependencies are available.")
