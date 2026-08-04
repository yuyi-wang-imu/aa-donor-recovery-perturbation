#!/usr/bin/env Rscript

# Submission-layout rebuild of Supplementary Figures S1-S11.
# This script only composes frozen, already-reviewed figure assets. It does not
# recompute statistics, change numerical results, or overwrite source files.

suppressPackageStartupMessages({
  library(magick)
  library(officer)
  library(ggplot2)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript build_supplementary_figures_typography_final_20260726_v2.R <ascii_input_dir> <ascii_output_dir>")
}

input_dir <- args[[1]]
if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

W <- 4320L
stamp <- "20260726_v4"

read_img <- function(name) {
  image_read(file.path(input_dir, name)) |>
    image_background("white", flatten = TRUE) |>
    image_convert(colorspace = "sRGB")
}

crop_px <- function(img, x, y, w, h) {
  info <- image_info(img)[1, ]
  x <- max(0L, min(as.integer(x), as.integer(info$width) - 1L))
  y <- max(0L, min(as.integer(y), as.integer(info$height) - 1L))
  w <- max(1L, min(as.integer(w), as.integer(info$width) - x))
  h <- max(1L, min(as.integer(h), as.integer(info$height) - y))
  image_crop(img, sprintf("%dx%d+%d+%d", w, h, x, y), repage = TRUE)
}

fit_width <- function(img, width = W) image_resize(img, sprintf("%dx", width))

pad_height <- function(img, height, gravity = "center") {
  info <- image_info(img)[1, ]
  image_extent(img, sprintf("%dx%d", info$width, height),
               gravity = gravity, color = "white")
}

tile_h <- function(imgs, width = W, gap = 40L) {
  each_w <- floor((width - gap * (length(imgs) - 1L)) / length(imgs))
  imgs <- lapply(imgs, fit_width, width = each_w)
  heights <- vapply(imgs, function(x) image_info(x)$height[[1]], numeric(1))
  imgs <- lapply(imgs, pad_height, height = max(heights))
  spacer <- image_blank(gap, max(heights), "white")
  parts <- list()
  for (i in seq_along(imgs)) {
    parts[[length(parts) + 1L]] <- imgs[[i]]
    if (i < length(imgs)) parts[[length(parts) + 1L]] <- spacer
  }
  image_append(image_join(parts), stack = FALSE)
}

tile_v <- function(imgs, width = W, gap = 42L) {
  imgs <- lapply(imgs, fit_width, width = width)
  spacer <- image_blank(width, gap, "white")
  parts <- list()
  for (i in seq_along(imgs)) {
    parts[[length(parts) + 1L]] <- imgs[[i]]
    if (i < length(imgs)) parts[[length(parts) + 1L]] <- spacer
  }
  image_append(image_join(parts), stack = TRUE)
}

info_band <- function(title, lines, width = W, height = 590L) {
  band <- image_blank(width, height, "white")
  band <- image_annotate(
    band, title, location = "+90+68", gravity = "northwest",
    size = 60, font = "Arial", weight = 700, color = "#172433"
  )
  y <- 175
  for (line in lines) {
    band <- image_annotate(
      band, line, location = sprintf("+95+%d", y), gravity = "northwest",
      size = 42, font = "Arial", color = "#44566c"
    )
    y <- y + 78
  }
  band
}

add_label <- function(img, label) {
  image_annotate(
    img, label, location = "+38+28", gravity = "northwest",
    size = 72, font = "Arial", weight = 700, color = "#111111"
  )
}

render_s8_robustness <- function(path) {
  robustness <- data.frame(
    module = factor(
      c("black", "blue", "cyan", "magenta", "purple", "turquoise", "yellow",
        "pink", "red", "royalblue", "darkgreen"),
      levels = c("black", "blue", "cyan", "magenta", "purple", "turquoise",
                 "yellow", "pink", "red", "royalblue", "darkgreen")
    ),
    consistent = c(9L, 9L, 9L, 9L, 9L, 9L, 9L, 8L, 8L, 8L, 7L)
  )
  robustness$percent <- 100 * robustness$consistent / 9
  robustness$class <- factor(
    ifelse(robustness$consistent == 9L, "9 of 9 settings",
           ifelse(robustness$consistent == 8L, "8 of 9 settings", "7 of 9 settings")),
    levels = c("9 of 9 settings", "8 of 9 settings", "7 of 9 settings")
  )
  robustness$label <- paste(robustness$consistent, "of 9")

  p <- ggplot(robustness, aes(module, percent, fill = class)) +
    geom_col(width = 0.82) +
    geom_text(aes(label = label), vjust = -0.55, size = 7.3, family = "Arial") +
    geom_hline(yintercept = 80, linetype = "dashed", linewidth = 0.8,
               colour = "#8a8a8a") +
    scale_fill_manual(
      values = c("9 of 9 settings" = "#3F8248",
                 "8 of 9 settings" = "#F0AC35",
                 "7 of 9 settings" = "#D34B5D"),
      name = NULL
    ) +
    scale_y_continuous(
      limits = c(0, 106),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = "Robustness of bone-marrow module localization",
      x = NULL,
      y = "Agreement across nine analysis settings"
    ) +
    theme_classic(base_family = "Arial", base_size = 24) +
    theme(
      plot.title = element_text(face = "bold", size = 31, hjust = 0.5,
                                margin = margin(b = 24)),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 21),
      axis.text.y = element_text(size = 21),
      axis.title.y = element_text(size = 24, margin = margin(r = 18)),
      legend.position = "bottom",
      legend.text = element_text(size = 21),
      legend.key.width = grid::unit(1.0, "cm"),
      plot.margin = margin(34, 42, 28, 44)
    )

  ragg::agg_png(path, width = 3600, height = 1700, units = "px",
                res = 300, background = "white")
  print(p)
  dev.off()
  invisible(path)
}

# Frozen sources. The legacy supplementary panels are used only where their
# scientific content has not changed. Corrected current main-figure assets are
# used for all terminology, P-value, panel-label, gene-format and 100-ns fixes.
old_s1 <- read_img("legacy_s1.png")
old_s4 <- read_img("legacy_s4.png")
old_s5 <- read_img("legacy_s5.png")
old_s8 <- read_img("legacy_s8.png")
old_s10 <- read_img("legacy_s10.png")

f1 <- read_img("current_figure1.png")
f2 <- read_img("current_figure2.png")
f3 <- read_img("current_figure3.png")
f5 <- read_img("current_figure5.png")
cd34 <- read_img("current_cd34_context.png")
md100 <- read_img("current_md100_fullmatrix.png")
s11 <- read_img("current_s11.png")

figs <- vector("list", 11L)

# S1: frozen traceable retrieval/manual-curation workflow.
figs[[1]] <- fit_width(old_s1, W)

# S2: complete prescription-mining evidence, including all 12 frozen rules.
figs[[2]] <- fit_width(f1, W)

# S3: candidate expression/marker-program panels plus detection coverage.
f2_info <- image_info(f2)[1, ]
s3_top <- crop_px(f2, 0, 0, f2_info$width, floor(f2_info$height * 0.50))
cd_info <- image_info(cd34)[1, ]
s3_det <- crop_px(cd34, floor(cd_info$width * 0.50), 0,
                  ceiling(cd_info$width * 0.50), floor(cd_info$height * 0.50))
s3_det <- image_extent(fit_width(s3_det, 3450), sprintf("%dx%d", W, 1880),
                       gravity = "center", color = "white")
figs[[3]] <- tile_v(list(s3_top, s3_det), W)

# S4: preserve PCA/correlation, replace the lower bulk panel with the corrected
# CD34 superscript/gene-italic/P-value version derived from the same frozen data.
s4_info <- image_info(old_s4)[1, ]
s4_top <- crop_px(old_s4, 0, 0, s4_info$width, floor(s4_info$height * 0.38))
s4_bulk <- crop_px(cd34, floor(cd_info$width * 0.50), floor(cd_info$height * 0.50),
                   ceiling(cd_info$width * 0.50), ceiling(cd_info$height * 0.50))
s4_bulk <- image_extent(fit_width(s4_bulk, 3250), sprintf("%dx%d", W, 2250),
                        gravity = "center", color = "white")
figs[[4]] <- tile_v(list(s4_top, s4_bulk), W)

# S5: retain the validated sampling and soft-threshold panels; replace only the
# explanatory band to standardize CD34-positive terminology and mathematical typography.
s5_info_old <- image_info(old_s5)[1, ]
s5_top <- crop_px(old_s5, 0, 0, s5_info_old$width, floor(s5_info_old$height * 0.73))
s5_band <- info_band(
  "Verified CD34-positive HSPC network inputs and parameters",
  c(
    "48 sequencing-library profiles; six technical-replicate sets aggregated at the raw-count level; 42 participant-by-time-point profiles from 23 participants.",
    "Signed bicor WGCNA; fixed soft-threshold power = 18; scale-free topology fit R² = 0.8016; mean connectivity = 19.3545.",
    "Minimum module size = 30; merge cut height = 0.25; 27 modules including grey; participant blocking retained in the predefined contrasts."
  ),
  height = 560
)
figs[[5]] <- tile_v(list(s5_top, s5_band), W)

# S6: corrected module contrasts and formal candidate evidence annotation.
f2_info <- image_info(f2)[1, ]
s6_main <- crop_px(f2, 0, floor(f2_info$height * 0.47),
                   f2_info$width, ceiling(f2_info$height * 0.53))
s6_band <- info_band(
  "Network-partition sensitivity",
  c(
    "The adjusted Rand index comparing the 48-library and 42-participant-by-time-point module partitions was 0.3948.",
    "Module colours are algorithmic labels. Module-level support is not interpreted as gene-level differential expression."
  ),
  height = 430
)
figs[[6]] <- tile_v(list(s6_main, s6_band), W)

# S7: corrected atlas/composition panels and explicit merged-category notation.
f3_info <- image_info(f3)[1, ]
s7_main <- crop_px(f3, 0, 0, f3_info$width, floor(f3_info$height * 0.50))
s7_band <- info_band(
  "Bone-marrow atlas integration and annotation summary",
  c(
    "768,617 cells were analysed; deterministic downsampling was used for display only.",
    "Seven major compartments were retained: HSPC, erythroid, megakaryocyte, myeloid, T/NK, B/plasma and stromal/endothelial.",
    "A slash denotes a merged cell category; low-confidence cells were excluded from the compartment-fraction denominator."
  ),
  height = 560
)
figs[[7]] <- tile_v(list(s7_main, s7_band), W)

# S8: corrected 9-of-9 notation and italic gene symbols, plus a typography-only
# redraw of the frozen robustness counts from the same sensitivity analysis.
s8_main <- crop_px(f3, 0, floor(f3_info$height * 0.50),
                   f3_info$width, ceiling(f3_info$height * 0.50))
s8_render_path <- file.path(tempdir(), paste0("S8_robustness_", stamp, ".png"))
render_s8_robustness(s8_render_path)
s8_sens <- image_read(s8_render_path) |>
  image_background("white", flatten = TRUE) |>
  image_convert(colorspace = "sRGB")
s8_sens <- image_extent(fit_width(s8_sens, 3250), sprintf("%dx%d", W, 1750),
                        gravity = "center", color = "white")
figs[[8]] <- tile_v(list(s8_main, s8_sens), W)

# S9: enlarged network panels and explicit edge semantics.
s9_band <- info_band(
  "Node and edge definitions",
  c(
    "Left network: herbs, resource-reported compounds and the formal 30-gene candidate order.",
    "Right network: candidate genes, WGCNA modules, biological programmes, marrow compartments and clinical contrasts form separate evidence layers.",
    "Edges denote recorded associations or evidence mappings; edge presence does not establish causality or treatment efficacy."
  ),
  height = 590
)
figs[[9]] <- tile_v(list(fit_width(f5, W), s9_band), W)

# S10: preserve the complete docking matrix and representative three-dimensional
# poses, remove the legacy short-trajectory detail row, and append the standardized
# five-complex 100-ns descriptor matrix.
s10_info_old <- image_info(old_s10)[1, ]
s10_docking <- crop_px(old_s10, 0, 0, s10_info_old$width,
                       floor(s10_info_old$height * 0.65))
s10_band <- info_band(
  "Docking and molecular-dynamics interpretation boundary",
  c(
    "The docking matrix retains all 138 receptor-ligand Vina mode-1 scores; positive scores remain visible and are not treated as structural support.",
    "The standardized molecular-dynamics matrix uses complete 0-100 ns trajectories for all five representative complexes and identical metric definitions."
  ),
  height = 520
)
figs[[10]] <- tile_v(list(s10_docking, fit_width(md100, W), s10_band), W)

# S11: current matched-control perturbation sensitivity output.
figs[[11]] <- fit_width(s11, W)

captions <- c(
  "Supplementary Figure S1. Literature retrieval, manual assessment, prescription-level extraction and herb-name standardization for the aplastic-anemia prescription dataset. Database-specific Chinese and English search strategies were applied to CNKI, Wanfang, PubMed and Web of Science through August 2025. Records were managed in NoteExpress and assessed manually by the researchers. The final analytical unit was an independently identifiable prescription record; 390 prescription records formed the frozen dataset used for subsequent analyses.",
  "Supplementary Figure S2. Prescription co-occurrence matrix and network for the 30 most frequent herbs, between-groups-linkage clustering and the complete two-herb-antecedent association rules. The 12 rules are presented in their frozen confidence-ranked order; point area denotes joint support and colour denotes lift.",
  "Supplementary Figure S3. Detection fractions and expression distributions of the 126 candidate genes in CD34-positive hematopoietic stem and progenitor cell profiles. All 126 genes were mapped, 125 were detected, and group-level expression and marker-program relationships are shown across 42 participant-by-time-point profiles from 23 participants.",
  "Supplementary Figure S4. Quality control of the GSE165870 lineage-negative CD34-positive transcriptomic dataset and between-group expression of the 126 candidates. Healthy and aplastic-anemia samples are shown by principal-component analysis and sample correlation, together with the independent bulk-expression context for the prespecified candidate set.",
  "Supplementary Figure S5. Technical-replicate aggregation, gene filtering and fixed soft-threshold diagnostics for participant-by-time-point CD34-positive HSPC profiles. Six technical-replicate sets were aggregated at the raw-count level, yielding 42 participant-by-time-point profiles from 23 participants. The signed bicor WGCNA used fixed power 18.",
  "Supplementary Figure S6. Co-expression-module statistics across the three predefined contrasts, evidence annotation of the 30 priority candidates and sensitivity of network partitioning to technical-replicate aggregation. Module-level support is distinguished from gene-level differential expression.",
  "Supplementary Figure S7. Integration, quality control and major-compartment annotation of the bone-marrow single-cell dataset. The atlas includes 768,617 cells and seven retained major marrow compartments. T/NK, B/plasma and stromal/endothelial denote merged cell categories; deterministic downsampling was used only for display.",
  "Supplementary Figure S8. Sensitivity analysis of candidate-associated module projections across marrow compartments and expression of genes in the THPO–MPL signaling axis. Concordance is reported as the number of consistent results among nine analysis settings, using the forms 9 of 9, 8 of 9 and 7 of 9.",
  "Supplementary Figure S9. Complete node and edge definitions and enlarged details for the herb-compound-candidate and candidate-module-programme-compartment networks. Network edges encode recorded associations or evidence mappings and are not interpreted as causal effects.",
  "Supplementary Figure S10. Complete docking-score matrix, representative predicted poses and standardized 100-ns molecular-dynamics descriptors for five representative complexes. The docking matrix retains all 138 receptor-ligand calculations, including three positive Vina mode-1 scores that are excluded from structural-support interpretation. Molecular-dynamics panels use the complete 0–100 ns trajectories and identical metric definitions across complexes.",
  "Supplementary Figure S11. scTenifoldKnk perturbation results under matched-control and cross-scenario sensitivity analyses. The outputs represent model-derived network responses and do not constitute experimental knockout validation."
)

png_paths <- character(length(figs))
for (i in seq_along(figs)) {
  img <- figs[[i]] |>
    image_background("white", flatten = TRUE) |>
    image_convert(colorspace = "sRGB") |>
    fit_width(W)
  base <- sprintf("Supplementary_Figure_S%d_%s", i, stamp)
  png_paths[[i]] <- file.path(output_dir, paste0(base, ".png"))
  image_write(img, png_paths[[i]], format = "png", density = "600x600")
  image_write(img, file.path(output_dir, paste0(base, ".tiff")),
              format = "tiff", density = "600x600", compression = "lzw")
  image_write(img, file.path(output_dir, paste0(base, ".pdf")),
              format = "pdf", density = "600x600")
}

pdf_pages <- lapply(figs, function(x) {
  x |>
    image_background("white", flatten = TRUE) |>
    image_convert(colorspace = "sRGB") |>
    fit_width(W)
})
image_write(
  image_join(pdf_pages),
  file.path(output_dir, sprintf("Supplementary_Figures_S1-S11_%s.pdf", stamp)),
  format = "pdf", density = "600x600"
)

normal_fp <- fp_text(
  font.family = "Times New Roman", font.size = 10,
  color = "000000", bold = FALSE
)
bold_fp <- fp_text(
  font.family = "Times New Roman", font.size = 10,
  color = "000000", bold = TRUE
)
heading_fp <- fp_text(
  font.family = "Times New Roman", font.size = 14,
  color = "000000", bold = TRUE
)

add_caption <- function(doc, i, caption) {
  prefix <- sprintf("Supplementary Figure S%d.", i)
  remainder <- sub(paste0("^", prefix, "\\s*"), "", caption)
  body_add_fpar(
    doc,
    fpar(
      ftext(prefix, bold_fp),
      ftext(paste0(" ", remainder), normal_fp),
      fp_p = fp_par(text.align = "left", line_spacing = 1.0)
    )
  )
}

cap_doc <- read_docx()
cap_doc <- body_add_fpar(cap_doc, fpar(ftext("Supplementary Figure Legends", heading_fp)))
for (i in seq_along(captions)) cap_doc <- add_caption(cap_doc, i, captions[[i]])
print(
  cap_doc,
  target = file.path(output_dir, sprintf("Supplementary_Figure_Legends_S1-S11_%s.docx", stamp))
)

review_doc <- read_docx()
review_doc <- body_add_fpar(review_doc, fpar(ftext("Supplementary Figures S1-S11", heading_fp)))
for (i in seq_along(figs)) {
  inf <- image_info(image_read(png_paths[[i]]))[1, ]
  h_in <- min(8.25, 6.65 * inf$height / inf$width)
  review_doc <- body_add_img(review_doc, src = png_paths[[i]], width = 6.65, height = h_in)
  review_doc <- add_caption(review_doc, i, captions[[i]])
  if (i < length(figs)) review_doc <- body_add_break(review_doc)
}
print(
  review_doc,
  target = file.path(output_dir, sprintf("Supplementary_Figures_S1-S11_review_copy_%s.docx", stamp))
)

thumbs <- lapply(seq_along(figs), function(i) {
  x <- figs[[i]] |> fit_width(1000)
  x <- image_border(x, "#b8c2cc", "2x2")
  add_label(x, paste0("S", i))
})
rows <- lapply(split(thumbs, ceiling(seq_along(thumbs) / 3)), function(z) {
  tile_h(z, width = 3100, gap = 35)
})
contact <- tile_v(rows, width = 3100, gap = 35)
image_write(
  contact,
  file.path(output_dir, sprintf("Supplementary_Figures_S1-S11_contact_sheet_%s.png", stamp)),
  format = "png", density = "300x300"
)

cat("Built", length(figs), "supplementary figures in", output_dir, "\n")
