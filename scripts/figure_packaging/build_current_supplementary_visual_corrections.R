#!/usr/bin/env Rscript

# Visual-only, non-overwriting corrections for Additional file 2.
# Scientific values and all unaffected pages remain unchanged.
# Corrected pages:
#   S6: redraw the reviewed 9/9, 8/9 and 7/9 robustness counts with the
#       legend outside the plotting region.
#   S7: remove the redundant page-level title while retaining the complete
#       138-pair matrix title, values and annotations.

suppressPackageStartupMessages({
  library(ggplot2)
  library(magick)
  library(pdftools)
  library(qpdf)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript build_v19_supplementary_visual_corrections_R_20260821_v1.R <source_pdf> <output_dir>", call. = FALSE)
}

source_pdf <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
if (dir.exists(output_dir)) stop("Output directory already exists: ", output_dir, call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

page_sizes <- pdf_pagesize(source_pdf)
if (nrow(page_sizes) != 16L) stop("Expected 16 source pages", call. = FALSE)

render_page <- function(page, dpi = 300L) {
  image_read(pdf_render_page(source_pdf, page = page, dpi = dpi)) |>
    image_background("white", flatten = TRUE) |>
    image_convert(colorspace = "sRGB")
}

write_raster_pdf <- function(img, path, dpi = 300L) {
  image_write(
    img, path, format = "pdf", density = sprintf("%dx%d", dpi, dpi),
    quality = 96, compression = "jpeg"
  )
  invisible(path)
}

# ---- S6: reviewed robustness counts; move legend fully above the bars. ----
robustness <- data.frame(
  module = factor(
    c("black", "blue", "cyan", "magenta", "purple", "turquoise", "yellow", "pink", "red", "royalblue", "darkgreen"),
    levels = c("black", "blue", "cyan", "magenta", "purple", "turquoise", "yellow", "pink", "red", "royalblue", "darkgreen")
  ),
  consistent = c(9L, 9L, 9L, 9L, 9L, 9L, 9L, 8L, 8L, 8L, 7L)
)
robustness$percent <- 100 * robustness$consistent / 9
robustness$class <- factor(
  ifelse(robustness$consistent == 9L, "100% agreement",
         ifelse(robustness$consistent == 8L, ">=80% agreement", "<80% agreement")),
  levels = c("100% agreement", ">=80% agreement", "<80% agreement")
)
robustness$label <- paste0(robustness$consistent, "/9")

p6 <- ggplot(robustness, aes(module, percent, fill = class)) +
  geom_col(width = 0.82) +
  geom_text(aes(label = label), vjust = -0.45, size = 3.45, family = "Arial") +
  geom_hline(yintercept = 80, linetype = "dashed", linewidth = 0.45, colour = "#8A8A8A") +
  scale_fill_manual(
    values = c("100% agreement" = "#3F8248", ">=80% agreement" = "#F0AC35", "<80% agreement" = "#D34B5D"),
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(0, 108), breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Supplementary Figure S6",
    subtitle = "Robustness of bone-marrow module localization",
    x = NULL,
    y = "Top-cell-type agreement across sensitivity variants (%)"
  ) +
  theme_classic(base_family = "Arial", base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 14.5, hjust = 0, margin = margin(b = 3)),
    plot.subtitle = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 5)),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 9.5, face = "bold"),
    axis.text.y = element_text(size = 9.5),
    axis.title.y = element_text(size = 10.2, margin = margin(r = 8)),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = 9.2),
    legend.key.width = grid::unit(0.55, "cm"),
    legend.margin = margin(t = 0, r = 0, b = 4, l = 0),
    plot.margin = margin(12, 20, 12, 18)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

s6_pdf <- file.path(output_dir, "page_06_S6_corrected.pdf")
ggsave(s6_pdf, p6, device = cairo_pdf, width = 8, height = 5.125, units = "in")

# ---- S7: retain all matrix pixels, remove only the redundant top-right title. ----
s7 <- render_page(7L, 300L)
s7_info <- image_info(s7)[1, ]
s7_draw <- image_draw(s7)
par(mar = c(0, 0, 0, 0), usr = c(0, s7_info$width, s7_info$height, 0), xpd = NA)
rect(s7_info$width * 0.80, 15, s7_info$width - 5, 100, col = "white", border = NA)
dev.off()
s7_pdf <- file.path(output_dir, "page_07_S7_corrected.pdf")
write_raster_pdf(s7_draw, s7_pdf, 300L)

# ---- Assemble without altering pages 1-5 or 8-16. ----
part_1_5 <- file.path(output_dir, "pages_01_05.pdf")
part_8_16 <- file.path(output_dir, "pages_08_16.pdf")
pdf_subset(source_pdf, pages = 1:5, output = part_1_5)
pdf_subset(source_pdf, pages = 8:16, output = part_8_16)

final_pdf <- file.path(output_dir, "Additional_file_2_Supplementary_Figures_S1-S16_v19_visual_corrected.pdf")
pdf_combine(c(part_1_5, s6_pdf, s7_pdf, part_8_16), output = final_pdf)

final_info <- pdf_info(final_pdf)
if (final_info$pages != 16L) stop("Final PDF does not contain 16 pages", call. = FALSE)
if (file.info(final_pdf)$size >= 20000000) stop("Final PDF is not below the 20,000,000-byte gate", call. = FALSE)

# ---- Pixel-level regression for all unaffected pages at 180 dpi. ----
render_dir <- file.path(output_dir, "pages_180dpi")
dir.create(render_dir, showWarnings = FALSE)
unaffected <- setdiff(1:16, c(6L, 7L))
qa <- data.frame(page = 1:16, status = "corrected_and_review_required", stringsAsFactors = FALSE)
for (i in 1:16) {
  before <- pdf_render_page(source_pdf, page = i, dpi = 180)
  after <- pdf_render_page(final_pdf, page = i, dpi = 180)
  image_write(image_read(after), file.path(render_dir, sprintf("page_%02d.png", i)), format = "png")
  if (i %in% unaffected) {
    identical_pixels <- identical(as.integer(before), as.integer(after))
    qa$status[[i]] <- if (identical_pixels) "unchanged_pixel_identical" else "FAIL_unexpected_change"
  }
}
if (any(grepl("^FAIL", qa$status))) stop("At least one unaffected page changed", call. = FALSE)
write.table(qa, file.path(output_dir, "page_regression_QA.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

sizes <- pdf_pagesize(final_pdf)
write.table(sizes, file.path(output_dir, "page_sizes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("FINAL_PDF=", final_pdf, "\n", sep = "")
cat("BYTES=", file.info(final_pdf)$size, "\n", sep = "")
cat("PAGES=", final_info$pages, "\n", sep = "")
