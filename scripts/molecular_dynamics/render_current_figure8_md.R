#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ragg)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript scripts/molecular_dynamics/render_current_figure8_md.R <new-output-directory>")
}

root <- getwd()
if (!file.exists(file.path(root, "CURRENT_MANUSCRIPT_ASSET_CHECKSUMS.tsv"))) {
  stop("Run this command from the repository root")
}
out_dir <- args[[1L]]
if (dir.exists(out_dir) || file.exists(out_dir)) stop("Refusing to overwrite existing output: ", out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source_dir <- file.path("derived_data", "molecular_dynamics", "current_figure8")
time_path <- file.path(source_dir, "Figure8_five_candidates_time_series_source.tsv.gz")
rmsf_path <- file.path(source_dir, "Figure8_five_candidates_ca_rmsf_source.tsv")
summary_path <- file.path(source_dir, "Figure8_final20ns_quantitative_summary.tsv")
inputs <- c(time_path, rmsf_path, summary_path)
if (any(!file.exists(inputs))) stop("Missing current Figure 8 source: ", paste(inputs[!file.exists(inputs)], collapse = "; "))
if (any(file.info(inputs)$size <= 0)) stop("Empty current Figure 8 source detected")

read_tsv <- function(path) {
  if (grepl("[.]gz$", path)) {
    connection <- gzfile(path, open = "rt")
    on.exit(close(connection), add = TRUE)
    read.delim(connection, check.names = FALSE)
  } else {
    read.delim(path, check.names = FALSE)
  }
}

ts <- read_tsv(time_path)
rmsf <- read_tsv(rmsf_path)
summary <- read_tsv(summary_path)

case_order <- c(
  "TOP2A_sesamin",
  "GSK3B_linarin",
  "KIT_3O_Methylorobol",
  "HIF1A_ARNT_butin_exploratory",
  "SYK_isofucosterol"
)
display <- c(
  TOP2A_sesamin = "TOP2A\u2013sesamin",
  GSK3B_linarin = "GSK3B\u2013linarin",
  KIT_3O_Methylorobol = "KIT\u20133\u2032-O-methylorobol",
  HIF1A_ARNT_butin_exploratory = "HIF1A\u2013ARNT\u2013butin\u2020",
  SYK_isofucosterol = "SYK\u2013isofucosterol"
)
colors <- c(
  TOP2A_sesamin = "#00897B",
  GSK3B_linarin = "#D97706",
  KIT_3O_Methylorobol = "#2563A6",
  HIF1A_ARNT_butin_exploratory = "#7C5AA6",
  SYK_isofucosterol = "#607D8B"
)
metric_order <- c(
  "protein_backbone_rmsd",
  "ligand_heavy_atom_rmsd_after_protein_fit",
  "global_protein_ligand_minimum_distance",
  "ligand_heavy_atom_proximity_fraction_lt_0p45nm"
)

required_ts <- c("case_id", "series", "time_ns", "value", "smooth")
required_rmsf <- c("case_id", "residue_index", "chain", "rmsf_nm")
required_summary <- c(
  "case_id", "protein_rmsd_mean_nm", "ligand_rmsd_mean_nm",
  "minimum_distance_median_nm", "proximity_fraction_mean", "ca_rmsf_median_nm"
)
if (!all(required_ts %in% names(ts))) stop("Time-series source columns are incomplete")
if (!all(required_rmsf %in% names(rmsf))) stop("RMSF source columns are incomplete")
if (!all(required_summary %in% names(summary))) stop("Final-20-ns summary columns are incomplete")
if (!identical(sort(unique(ts$case_id)), sort(case_order))) stop("Unexpected time-series case membership")
if (!identical(sort(unique(rmsf$case_id)), sort(case_order))) stop("Unexpected RMSF case membership")
if (!identical(sort(unique(summary$case_id)), sort(case_order))) stop("Unexpected summary case membership")
if (any(grepl("PARP1|CDK6", c(ts$case_id, rmsf$case_id, summary$case_id)))) stop("Excluded legacy cases detected")
if (any(!is.finite(ts$time_ns)) || any(!is.finite(ts$value))) stop("Non-finite time-series value")
if (any(!is.finite(rmsf$residue_index)) || any(!is.finite(rmsf$rmsf_nm)) || any(rmsf$rmsf_nm < 0)) stop("Invalid RMSF value")

for (case in case_order) {
  for (metric in metric_order) {
    block <- ts[ts$case_id == case & ts$series == metric, ]
    if (nrow(block) != 10001L) stop(case, "/", metric, ": expected 10001 points")
    block <- block[order(block$time_ns), ]
    if (abs(block$time_ns[1L]) > 0.02 || abs(tail(block$time_ns, 1L) - 100) > 0.02) stop(case, "/", metric, ": wrong time span")
    if (any(diff(block$time_ns) <= 0) || abs(median(diff(block$time_ns)) - 0.01) > 1e-4) stop(case, "/", metric, ": nonuniform spacing")
  }
}
proximity <- ts$value[ts$series == metric_order[[4L]]]
if (any(proximity < 0 | proximity > 1)) stop("Proximity fraction outside [0, 1]")

ts$case_id <- factor(ts$case_id, levels = case_order)
rmsf$case_id <- factor(rmsf$case_id, levels = case_order)
summary$case_id <- factor(summary$case_id, levels = case_order)
ts <- ts[order(ts$case_id, ts$series, ts$time_ns), ]
rmsf <- rmsf[order(rmsf$case_id, rmsf$residue_index), ]
summary <- summary[match(case_order, as.character(summary$case_id)), ]

base_theme <- theme_classic(base_family = "Arial", base_size = 7.2) +
  theme(
    plot.title = element_text(size = 8.4, face = "bold", color = "black", hjust = 0),
    axis.title = element_text(size = 7.2, color = "black"),
    axis.text = element_text(size = 6.6, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    panel.grid.major = element_line(linewidth = 0.22, color = "#DDE3E8"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 6.4, color = "black"),
    legend.key.width = grid::unit(6.5, "mm"),
    plot.margin = margin(2.5, 3.5, 2.5, 3.5)
  )

metric_meta <- list(
  protein_backbone_rmsd = list(title = "Protein-backbone RMSD", y = "RMSD (nm)"),
  ligand_heavy_atom_rmsd_after_protein_fit = list(title = "Ligand RMSD after protein fit", y = "RMSD (nm)"),
  global_protein_ligand_minimum_distance = list(title = "Minimum protein\u2013ligand distance", y = "Distance (nm)"),
  ligand_heavy_atom_proximity_fraction_lt_0p45nm = list(title = "Ligand proximity fraction (<0.45 nm)", y = "Fraction")
)

make_ts_plot <- function(metric, tag) {
  d <- ts[ts$series == metric, ]
  p <- ggplot(d, aes(x = time_ns, y = value, group = case_id, color = case_id)) +
    annotate("rect", xmin = 80, xmax = 100, ymin = -Inf, ymax = Inf, fill = "#E9EDF1", alpha = 0.65) +
    geom_line(linewidth = 0.16, alpha = 0.10, na.rm = TRUE) +
    geom_line(aes(y = smooth), linewidth = 0.62, alpha = 0.98, na.rm = TRUE) +
    scale_color_manual(values = colors, breaks = case_order, labels = unname(display[case_order])) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = expansion(mult = c(0, 0))) +
    labs(title = paste0(tag, "   ", metric_meta[[metric]]$title), x = "Time (ns)", y = metric_meta[[metric]]$y) +
    base_theme
  if (metric == metric_order[[4L]]) {
    p <- p + scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0)))
  }
  p
}

pA <- make_ts_plot(metric_order[[1L]], "A")
pB <- make_ts_plot(metric_order[[2L]], "B") + theme(legend.position = "none")
pC <- make_ts_plot(metric_order[[3L]], "C") + theme(legend.position = "none")
pD <- make_ts_plot(metric_order[[4L]], "D") + theme(legend.position = "none")
top_grid <- (pA + pB) / (pC + pD) + plot_layout(guides = "collect") & theme(legend.position = "top")

make_rmsf <- function(case) {
  d <- rmsf[as.character(rmsf$case_id) == case, ]
  title <- unname(display[[case]])
  if (case == "TOP2A_sesamin") title <- paste0("E   ", title)
  p <- ggplot(d, aes(x = residue_index, y = rmsf_nm)) +
    geom_line(linewidth = 0.40, color = colors[[case]]) +
    labs(title = title, x = "Residue index", y = "C\u03b1 RMSF (nm)") +
    base_theme +
    theme(
      plot.title = element_text(size = 6.5, face = if (case %in% c("TOP2A_sesamin", "GSK3B_linarin")) "bold" else "plain", color = "black", hjust = 0.5),
      axis.title = element_text(size = 6.0),
      axis.text = element_text(size = 5.6),
      panel.grid.major = element_line(linewidth = 0.18, color = "#E4E8EB"),
      legend.position = "none",
      plot.margin = margin(2, 2, 2, 2)
    )
  if (case == "HIF1A_ARNT_butin_exploratory") {
    break_at <- max(d$residue_index[d$chain == "A"]) + 0.5
    p <- p + geom_vline(xintercept = break_at, linetype = "dashed", linewidth = 0.28, color = "#71717A")
  }
  p
}

rmsf_panel <- wrap_plots(lapply(case_order, make_rmsf), ncol = 5) + plot_annotation(tag_levels = NULL)

heat_long <- rbind(
  data.frame(case_id = as.character(summary$case_id), metric = "Protein RMSD\nmean (nm)", value = summary$protein_rmsd_mean_nm, digits = sprintf("%.3f", summary$protein_rmsd_mean_nm)),
  data.frame(case_id = as.character(summary$case_id), metric = "Ligand RMSD\nmean (nm)", value = summary$ligand_rmsd_mean_nm, digits = sprintf("%.3f", summary$ligand_rmsd_mean_nm)),
  data.frame(case_id = as.character(summary$case_id), metric = "Minimum distance\nmedian (nm)", value = summary$minimum_distance_median_nm, digits = sprintf("%.3f", summary$minimum_distance_median_nm)),
  data.frame(case_id = as.character(summary$case_id), metric = "Proximity fraction\nmean", value = summary$proximity_fraction_mean, digits = sprintf("%.3f", summary$proximity_fraction_mean))
)
heat_long$case_id <- factor(heat_long$case_id, levels = rev(case_order))
heat_long$display <- unname(display[as.character(heat_long$case_id)])
heat_levels <- unname(display[rev(case_order)])
heat_long$metric <- factor(heat_long$metric, levels = c("Protein RMSD\nmean (nm)", "Ligand RMSD\nmean (nm)", "Minimum distance\nmedian (nm)", "Proximity fraction\nmean"))
heat_long$scaled <- ave(heat_long$value, heat_long$metric, FUN = function(v) if (diff(range(v)) == 0) rep(0.5, length(v)) else (v - min(v)) / diff(range(v)))
heat_long$lead <- as.character(heat_long$case_id) %in% c("TOP2A_sesamin", "GSK3B_linarin")

pF <- ggplot(heat_long, aes(x = metric, y = factor(display, levels = heat_levels), fill = scaled)) +
  geom_tile(color = "white", linewidth = 0.75) +
  geom_text(aes(label = digits, fontface = ifelse(lead, "bold", "plain")), size = 2.45, color = "black") +
  scale_fill_gradient(low = "#F7FAFC", high = "#6E9EB8", limits = c(0, 1), guide = "none") +
  scale_x_discrete(position = "bottom") +
  labs(
    title = "F   Final 20 ns quantitative summary",
    x = NULL,
    y = NULL,
    caption = "Bold values: strongest donor-level and matched-background Geneformer support. \u2020 Exploratory HIF1A\u2013ARNT starting complex."
  ) +
  theme_minimal(base_family = "Arial", base_size = 7) +
  theme(
    plot.title = element_text(size = 8.4, face = "bold", color = "black", hjust = 0),
    axis.text.x = element_text(size = 6.4, color = "black", lineheight = 0.95),
    axis.text.y = element_text(size = 6.5, color = "black"),
    panel.grid = element_blank(),
    plot.caption = element_text(size = 5.8, color = "#343A40", hjust = 0),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 4, 2, 4)
  )

final_plot <- wrap_plots(
  list(top_grid, rmsf_panel, pF),
  ncol = 1,
  heights = c(3.25, 1.15, 1.25)
) & theme(plot.background = element_rect(fill = "white", color = NA))

png_path <- file.path(out_dir, "Figure_8.png")
ggsave(
  png_path,
  final_plot,
  device = ragg::agg_png,
  width = 183,
  height = 225,
  units = "mm",
  dpi = 600,
  bg = "white"
)

qa <- data.frame(
  check = c(
    "five_current_cases",
    "legacy_cases_absent",
    "10001_frames_per_case_and_metric",
    "zero_to_100_ns",
    "proximity_in_unit_interval",
    "output_nonempty"
  ),
  result = c(
    "PASS",
    "PASS",
    if (all(table(ts$case_id, ts$series) == 10001L)) "PASS" else "FAIL",
    "PASS",
    "PASS",
    if (file.exists(png_path) && file.info(png_path)$size > 0) "PASS" else "FAIL"
  ),
  stringsAsFactors = FALSE
)
write.table(qa, file.path(out_dir, "Figure_8_render_QA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (any(qa$result != "PASS")) stop("Figure 8 render QA failed")

cat("PASS: current Figure 8 rendered to ", png_path, "\n", sep = "")
