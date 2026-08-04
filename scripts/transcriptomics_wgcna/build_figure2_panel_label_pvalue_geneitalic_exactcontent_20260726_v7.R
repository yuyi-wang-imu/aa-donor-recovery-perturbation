#!/usr/bin/env Rscript

# Figure contract
# Core conclusion: candidate-gene expression and marker-program correlations
# provide CD34-positive HSPC expression annotations; panels C-D retain the
# frozen analytical results already shown in the manuscript.
# Evidence boundary: this script changes typography, panel-label alignment,
# and terminology only. It does not alter values, clustering, candidate order,
# module evidence, or statistics. Backend: R only.

suppressPackageStartupMessages({
  library(grid)
  library(magick)
  library(pheatmap)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "Usage: Rscript script.R GROUP_SUMMARY.tsv CORRELATION.csv ",
    "PRIORITY.csv CURRENT_FIGURE2.png OUTPUT_DIR"
  )
}

group_summary_path <- args[[1L]]
correlation_path <- args[[2L]]
priority_path <- args[[3L]]
current_figure_path <- args[[4L]]
output_dir <- args[[5L]]
for (path in c(group_summary_path, correlation_path, priority_path, current_figure_path)) {
  if (!file.exists(path)) stop("Required input not found: ", path)
}
if (!dir.exists(output_dir)) stop("Output directory not found: ", output_dir)

priority_genes <- c("SYK", "TNF", "GSK3B", "AKT1", "PARP1")

scale_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

draw_pheatmap_with_italic_gene_labels <- function(...) {
  ph <- pheatmap(..., silent = TRUE)
  row_index <- which(ph$gtable$layout$name == "row_names")
  if (length(row_index) != 1L) {
    stop("Expected exactly one pheatmap row_names grob")
  }
  ph$gtable$grobs[[row_index]]$gp$fontface <- "italic"
  grid.newpage()
  grid.draw(ph$gtable)
}

group_summary <- read.delim(
  group_summary_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
priority <- read.csv(
  priority_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
target_panel_cor <- read.csv(
  correlation_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

group_expr_wide <- reshape(
  group_summary[, c("gene", "group", "mean_log1p_cp10k")],
  idvar = "gene",
  timevar = "group",
  direction = "wide"
)
names(group_expr_wide) <- sub("^mean_log1p_cp10k\\.", "", names(group_expr_wide))

top_expr_genes <- head(priority$gene, 55L)
heatmap_genes <- unique(c(
  priority_genes[priority_genes %in% group_expr_wide$gene],
  top_expr_genes
))
group_mat <- as.matrix(
  group_expr_wide[
    match(heatmap_genes, group_expr_wide$gene),
    c("HD", "SAA_baseline", "SAA_3M", "SAA_6M")
  ]
)
rownames(group_mat) <- heatmap_genes
group_mat[!is.finite(group_mat)] <- 0

top_cor_genes <- unique(c(
  priority_genes[priority_genes %in% target_panel_cor$gene],
  head(unique(target_panel_cor$gene), 45L)
))
cor_mat <- reshape(
  target_panel_cor[
    target_panel_cor$gene %in% top_cor_genes,
    c("gene", "panel", "spearman_rho")
  ],
  idvar = "gene",
  timevar = "panel",
  direction = "wide"
)
rownames(cor_mat) <- cor_mat$gene
cor_mat$gene <- NULL
names(cor_mat) <- sub("^spearman_rho\\.", "", names(cor_mat))
cor_mat <- as.matrix(cor_mat)
cor_mat[!is.finite(cor_mat)] <- 0

panel_a_path <- file.path(output_dir, "Figure2A_geneitalic_panel_label_v7_20260726.png")
panel_b_path <- file.path(output_dir, "Figure2B_geneitalic_panel_label_v7_20260726.png")
subtitle_c_path <- file.path(output_dir, "Figure2C_subtitle_pvalue_v7_20260726.png")
combined_path <- file.path(output_dir, "Figure2_full_exactcontent_geneitalic_panel_label_pvalue_v7_20260726.png")
for (path in c(panel_a_path, panel_b_path, subtitle_c_path, combined_path)) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
}

ragg::agg_png(
  panel_a_path,
  width = 2160,
  height = 1950,
  units = "px",
  res = 220,
  background = "white"
)
draw_pheatmap_with_italic_gene_labels(
  scale_rows(group_mat),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize = 13.5,
  fontsize_row = 7.1,
  fontsize_col = 12.5,
  angle_col = 0,
  labels_col = c("HD", "SAA baseline", "SAA 3M", "SAA 6M"),
  color = colorRampPalette(c("#26547C", "#F7F4EA", "#D1495B"))(101),
  main = "126 candidate genes in GSE247531 CD34\u207a HSPCs: group-level expression",
  border_color = "#8B949E"
)
dev.off()

panel_labels <- c(
  HSPC_identity = "HSPC identity",
  liver_THPO_MPL_response = "liver THPO-MPL response",
  interferon_inflammation = "interferon/inflammation",
  HSPC_injury_stress = "HSPC injury/stress",
  hematopoietic_support_context = "hematopoietic support"
)
display_labels <- unname(panel_labels[colnames(cor_mat)])
display_labels[is.na(display_labels)] <- gsub("_", " ", colnames(cor_mat)[is.na(display_labels)])

ragg::agg_png(
  panel_b_path,
  width = 2160,
  height = 1950,
  units = "px",
  res = 220,
  background = "white"
)
draw_pheatmap_with_italic_gene_labels(
  cor_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 13.5,
  fontsize_row = 7.1,
  fontsize_col = 11.0,
  angle_col = 45,
  labels_col = display_labels,
  color = colorRampPalette(c("#26547C", "#F7F4EA", "#D1495B"))(101),
  breaks = seq(-1, 1, length.out = 102),
  main = "Spearman correlation: candidate genes vs CD34 marker-program scores",
  border_color = "#8B949E"
)
dev.off()

# Render the only italic character in "P < 0.05" as a real italic glyph.
ragg::agg_png(
  subtitle_c_path,
  width = 1900,
  height = 120,
  units = "px",
  res = 220,
  background = "transparent"
)
par(mar = c(0, 0, 0, 0), family = "Arial")
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(
  0.002,
  0.50,
  bquote("Filled: BH FDR < 0.05; open: nominal " * italic(P) * " < 0.05"),
  adj = c(0, 0.5),
  cex = 1.15,
  col = "#5F6B76"
)
dev.off()

panel_box <- function(path, label) {
  img <- image_read(path)
  img <- image_trim(img, fuzz = 1)
  img <- image_resize(img, "2100x1950")
  img <- image_extent(img, "2160x2050", gravity = "northeast", color = "white")
  image_annotate(
    img,
    label,
    gravity = "northwest",
    location = "+10+24",
    size = 46,
    font = "Arial",
    weight = 700,
    color = "black"
  )
}

top <- image_append(
  c(panel_box(panel_a_path, "A"), panel_box(panel_b_path, "B")),
  stack = FALSE
)
top <- image_extent(top, "4320x2050", gravity = "center", color = "white")

current <- image_read(current_figure_path)
info <- image_info(current)
if (info$width != 4320L || info$height != 4300L) {
  stop("Unexpected current Figure 2 dimensions: ", info$width, "x", info$height)
}
bottom <- image_crop(current, "4320x2250+0+2050", repage = TRUE)

# Clear only the lower-panel title band, leaving all data marks unchanged.
title_band <- image_blank(width = 4320, height = 235, color = "white")
bottom <- image_composite(bottom, title_band, operator = "over", offset = "+0+0")
bottom <- image_annotate(
  bottom, "C", gravity = "northwest", location = "+10+24",
  size = 46, font = "Arial", weight = 700, color = "black"
)
bottom <- image_annotate(
  bottom, "Subject-blocked module eigengene contrasts",
  gravity = "northwest", location = "+100+10",
  size = 48, font = "Arial", weight = 700, color = "black"
)
bottom <- image_composite(
  bottom,
  image_read(subtitle_c_path),
  operator = "over",
  offset = "+100+92"
)
bottom <- image_annotate(
  bottom, "D", gravity = "northwest", location = "+2170+24",
  size = 46, font = "Arial", weight = 700, color = "black"
)
bottom <- image_annotate(
  bottom, "Predefined priority-candidate evidence annotation matrix",
  gravity = "northwest", location = "+2260+10",
  size = 48, font = "Arial", weight = 700, color = "black"
)
bottom <- image_annotate(
  bottom, "Formal 1–30 order retained; docking excluded",
  gravity = "northwest", location = "+2260+100",
  size = 34, font = "Arial", weight = 400, color = "#5F6B76"
)

combined <- image_append(c(top, bottom), stack = TRUE)
image_write(combined, combined_path, format = "png", density = "600x600")

cat(panel_a_path, "\n", panel_b_path, "\n", combined_path, "\n", sep = "")
