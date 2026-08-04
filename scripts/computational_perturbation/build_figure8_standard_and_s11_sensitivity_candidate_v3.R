#!/usr/bin/env Rscript

# Reproducible figures for two deliberately separated evidence layers:
# Figure 8 shows the standard ten-candidate scTenifoldKnk results, whereas
# Figure S11 shows the secondary 200-control matched sensitivity analysis.
# All panels are derived from the predefined 20260724 common-seed batch. The script
# refuses to overwrite any pre-existing output.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: Rscript build_figure8_standard_and_s11_sensitivity_candidate_v3.R ",
    "<derived_data_directory> <figure_output_directory>"
  )
}

data_dir <- normalizePath(args[[1]], mustWork = TRUE)
figure_dir <- normalizePath(args[[2]], mustWork = TRUE)
diff_dir <- file.path(
  data_dir, "current210_candidate_diffregulation_20260725"
)

paths <- list(
  manifest = file.path(data_dir, "run_manifest_210.csv"),
  endpoints = file.path(data_dir, "descriptive_endpoints_all_runs.csv"),
  assignments = file.path(
    data_dir, "global_unique_rematching_assignments_v1.csv"
  ),
  sensitivity = file.path(
    data_dir, "candidate_global_unique_calibration_sensitivity_v1.csv"
  ),
  original = file.path(data_dir, "candidate_matched_calibration_v3.csv"),
  pooled_null = file.path(
    data_dir, "pooled_global_unique_null_residuals_sensitivity_v1.csv"
  )
)
for (path in c(unlist(paths), diff_dir)) {
  if (!file.exists(path) && !dir.exists(path)) {
    stop("Missing source: ", path)
  }
}

core10 <- c(
  "CDK6", "CA2", "PARP1", "KIT", "SYK",
  "GSK3B", "HIF1A", "TOP2A", "TERT", "CD38"
)
target_colors <- c(
  CDK6 = "#2F6FA5", CA2 = "#4C8FC0", PARP1 = "#5D63A8",
  KIT = "#309B8C", SYK = "#7655A6", GSK3B = "#A8781D",
  HIF1A = "#8A5FA7", TOP2A = "#C2703D", TERT = "#8A929B",
  CD38 = "#3B7D78"
)
ink <- "#20242A"
muted <- "#65717D"
grid_col <- "#DDE3E8"
light_grey <- "#E8EDF1"
mid_grey <- "#82909C"
red <- "#B44949"
blue <- "#3F76A8"
teal <- "#33998B"
gold <- "#D99A18"

main_stem <- file.path(
  figure_dir,
  "Figure_8_standard_scTenifoldKnk_network_responses_candidate_20260725_v10"
)
supp_stem <- file.path(
  figure_dir,
  "Supplementary_Figure_S11_matched_control_sensitivity_candidate_20260725_v4"
)
main_outputs <- paste0(main_stem, c(".svg", ".pdf", ".png", ".tiff"))
supp_outputs <- paste0(supp_stem, c(".svg", ".pdf", ".png", ".tiff"))

source_paths <- c(
  file.path(
    data_dir,
    "Figure8_standard_panelB_response_landscape_source_20260725_v4.csv"
  ),
  file.path(
    data_dir,
    "Figure8_standard_panelC_ranked_response_source_20260725_v4.csv"
  ),
  file.path(
    data_dir,
    "Figure8_standard_panelD_response_breadth_source_20260725_v4.csv"
  ),
  file.path(
    data_dir,
    "FigureS11_matched_control_calibration_source_20260725_v4.csv"
  ),
  file.path(
    data_dir,
    "FigureS11_matched_control_covariate_balance_source_20260725_v4.csv"
  ),
  file.path(
    data_dir,
    "FigureS11_matched_control_pooled_null_source_20260725_v4.csv"
  )
)
if (any(file.exists(c(main_outputs, supp_outputs)))) {
  stop("Refusing to overwrite existing figure outputs")
}

manifest <- read.csv(paths$manifest, stringsAsFactors = FALSE)
endpoints <- read.csv(paths$endpoints, stringsAsFactors = FALSE)
assignments <- read.csv(paths$assignments, stringsAsFactors = FALSE)
sensitivity <- read.csv(paths$sensitivity, stringsAsFactors = FALSE)
original <- read.csv(paths$original, stringsAsFactors = FALSE)
pooled_null <- read.csv(paths$pooled_null, stringsAsFactors = FALSE)

if (nrow(assignments) != 200L || anyDuplicated(assignments$control_gene)) {
  stop("Expected 200 disjoint matched controls")
}
if (!identical(sensitivity$candidate, core10)) {
  stop("Candidate order mismatch in sensitivity results")
}
if (sum(sensitivity$individual_gate_pass) != 9L) {
  stop("Expected nine candidates to pass the strict matching gate")
}
if (any(sensitivity$BH_q < 0.05, na.rm = TRUE)) {
  stop("Unexpected calibrated BH q < 0.05")
}

diff_files <- list.files(
  diff_dir,
  pattern = "^[0-9]{3}_candidate_.*_diffRegulation[.]csv$",
  full.names = TRUE
)
if (length(diff_files) != 10L) {
  stop("Expected exactly 10 candidate differential-response files")
}
diff_list <- lapply(diff_files, function(path) {
  dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("gene", "distance", "Z", "FC", "p.value", "p.adj")
  if (!identical(names(dat), required) || nrow(dat) != 3000L) {
    stop("Unexpected candidate file structure: ", basename(path))
  }
  candidate <- sub(
    "^[0-9]{3}_candidate_(.*)_diffRegulation[.]csv$",
    "\\1", basename(path)
  )
  dat$candidate <- candidate
  dat
})
diff_all <- bind_rows(diff_list)
if (!setequal(unique(diff_all$candidate), core10)) {
  stop("Candidate symbols in differential-response files are incomplete")
}

diff_downstream <- diff_all |>
  filter(gene != candidate) |>
  mutate(
    significant = !is.na(p.adj) & p.adj < 0.05,
    neg_log10_q = pmin(-log10(pmax(p.adj, 1e-12)), 12)
  )

# A deterministic display subset: twelve genes with the largest significant
# absolute Z and twelve with the largest across-candidate Z variability.
gene_summary <- diff_downstream |>
  group_by(gene) |>
  summarise(
    significant_any = any(significant),
    max_abs_z_significant = ifelse(
      significant_any, max(abs(Z[significant]), na.rm = TRUE), NA_real_
    ),
    sd_z = sd(Z, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(significant_any)
top_peak <- gene_summary |>
  arrange(desc(max_abs_z_significant), gene) |>
  slice_head(n = 12) |>
  pull(gene)
top_variable <- gene_summary |>
  filter(!gene %in% top_peak) |>
  arrange(desc(sd_z), gene) |>
  slice_head(n = 12) |>
  pull(gene)
display_genes <- c(top_peak, top_variable)
if (length(display_genes) < 20L) {
  stop("Insufficient response genes for the display landscape")
}

heatmap_dat <- diff_downstream |>
  filter(gene %in% display_genes) |>
  select(candidate, gene, Z, p.adj, significant)

zmat <- heatmap_dat |>
  select(candidate, gene, Z) |>
  pivot_wider(names_from = candidate, values_from = Z) |>
  as.data.frame()
rownames(zmat) <- zmat$gene
zmat$gene <- NULL
zmat <- as.matrix(zmat[, core10, drop = FALSE])
row_order <- rownames(zmat)[hclust(dist(zmat), method = "average")$order]
heatmap_dat$candidate <- factor(heatmap_dat$candidate, levels = core10)
heatmap_dat$gene <- factor(heatmap_dat$gene, levels = rev(row_order))

rank_dat <- diff_downstream |>
  group_by(candidate) |>
  arrange(p.adj, desc(abs(Z)), .by_group = TRUE) |>
  mutate(response_rank = row_number()) |>
  ungroup() |>
  filter(response_rank <= 120) |>
  mutate(candidate = factor(candidate, levels = core10))

candidate_endpoints <- sensitivity |>
  transmute(
    candidate,
    candidate_endpoint,
    matched_control_median,
    candidate_standardized_residual,
    matched_control_local_percentile,
    max_abs_standardized_matching_imbalance,
    individual_gate_pass,
    empirical_calibration_p,
    BH_q,
    reporting_status
  )
candidate_endpoints$candidate <- factor(
  candidate_endpoints$candidate, levels = rev(core10)
)
if (any(candidate_endpoints$candidate_endpoint < 1L)) {
  stop("Every prespecified candidate must have at least one BH-significant response gene")
}

manifest_endpoints <- merge(
  manifest,
  endpoints[, c("run_id", "n_sig_excluding_gKO_padj_0_05")],
  by = "run_id", all.x = TRUE, sort = FALSE
)
control_rows <- lapply(core10, function(candidate_symbol) {
  genes <- assignments$control_gene[
    assignments$candidate == candidate_symbol
  ]
  rows <- manifest_endpoints[
    manifest_endpoints$run_role == "matched_control" &
      manifest_endpoints$gKO %in% genes,
    , drop = FALSE
  ]
  if (nrow(rows) != 20L) {
    stop("Control count mismatch for ", candidate_symbol)
  }
  data.frame(
    candidate = candidate_symbol,
    control_gene = rows$gKO,
    response_count = rows$n_sig_excluding_gKO_padj_0_05,
    stringsAsFactors = FALSE
  )
})
control_plot <- bind_rows(control_rows)
control_plot$candidate <- factor(control_plot$candidate, levels = rev(core10))

metric_map <- c(
  expression = "mean_log1p_cp10k",
  detection = "detected_fraction",
  out_degree = "out_degree",
  out_strength = "out_strength"
)
metric_labels <- c(
  expression = "Expression",
  detection = "Detection",
  out_degree = "Out-degree",
  out_strength = "Out-strength"
)
balance_rows <- list()
for (candidate_symbol in core10) {
  candidate_row <- manifest[
    manifest$run_role == "candidate" &
      manifest$gKO == candidate_symbol,
    , drop = FALSE
  ]
  control_genes <- assignments$control_gene[
    assignments$candidate == candidate_symbol
  ]
  controls <- manifest[
    manifest$run_role == "matched_control" &
      manifest$gKO %in% control_genes,
    , drop = FALSE
  ]
  if (nrow(candidate_row) != 1L || nrow(controls) != 20L) {
    stop("Manifest mismatch for ", candidate_symbol)
  }
  for (metric in names(metric_map)) {
    suffix <- metric_map[[metric]]
    candidate_value <- candidate_row[[paste0("candidate_", suffix)]][1]
    control_values <- controls[[paste0("control_", suffix)]]
    standardized_difference <- (
      candidate_value - mean(control_values)
    ) / sd(control_values)
    balance_rows[[length(balance_rows) + 1L]] <- data.frame(
      candidate = candidate_symbol,
      metric = metric_labels[[metric]],
      standardized_difference = standardized_difference,
      stringsAsFactors = FALSE
    )
  }
}
balance_metric <- bind_rows(balance_rows)
balance_metric$candidate <- factor(balance_metric$candidate, levels = rev(core10))
balance_metric$metric <- factor(
  balance_metric$metric,
  levels = c("Expression", "Detection", "Out-degree", "Out-strength")
)
check_max <- balance_metric |>
  mutate(candidate = as.character(candidate)) |>
  group_by(candidate) |>
  summarise(
    max_abs = max(abs(standardized_difference)),
    .groups = "drop"
  ) |>
  arrange(match(candidate, core10))
if (max(abs(
  check_max$max_abs -
    sensitivity$max_abs_standardized_matching_imbalance
)) > 1e-6) {
  stop("Recomputed covariate imbalance does not match the predefined source table")
}

write_or_verify_csv <- function(object, path, tolerance = 1e-8) {
  generated <- as.data.frame(
    object, stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!file.exists(path)) {
    write.csv(
      generated, path, row.names = FALSE, fileEncoding = "UTF-8"
    )
    cat("WROTE_SOURCE_DATA:", path, "\n")
    return(invisible(TRUE))
  }

  existing <- read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!identical(names(existing), names(generated)) ||
      nrow(existing) != nrow(generated)) {
    stop("Packaged source-data structure mismatch: ", path)
  }
  for (column in names(generated)) {
    expected <- generated[[column]]
    observed <- existing[[column]]
    if (is.numeric(expected) || is.integer(expected)) {
      same <- isTRUE(all.equal(
        as.numeric(observed), as.numeric(expected),
        tolerance = tolerance, check.attributes = FALSE
      ))
    } else {
      same <- identical(as.character(observed), as.character(expected))
    }
    if (!same) {
      stop("Packaged source-data value mismatch in ", column, ": ", path)
    }
  }
  cat("VERIFIED_EXISTING_SOURCE_DATA:", path, "\n")
  invisible(TRUE)
}

write_or_verify_csv(heatmap_dat, source_paths[[1]])
write_or_verify_csv(rank_dat, source_paths[[2]])
write_or_verify_csv(candidate_endpoints, source_paths[[3]])
write_or_verify_csv(sensitivity, source_paths[[4]])
write_or_verify_csv(balance_metric, source_paths[[5]])
write_or_verify_csv(pooled_null, source_paths[[6]])

theme_pub <- theme_minimal(base_size = 7.3, base_family = "sans") +
  theme(
    plot.title = element_text(
      face = "bold", size = 8.6, color = ink,
      margin = margin(b = 2.5)
    ),
    plot.subtitle = element_text(
      size = 6.4, color = muted, margin = margin(b = 4)
    ),
    axis.title = element_text(size = 6.9, color = ink),
    axis.text = element_text(size = 6.4, color = ink),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(
      color = grid_col, linewidth = 0.28
    ),
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6.1),
    plot.margin = margin(4, 5, 4, 5)
  )

panel_a <- ggplot() +
  annotate(
    "label", x = 1.05, y = 1,
    label = "10 prespecified\ncandidate genes",
    size = 2.45, linewidth = 0.3,
    color = ink, fill = "#E7EFF6"
  ) +
  annotate(
    "label", x = 3.35, y = 1,
    label = "CD34\u207a HSPC\nreference network",
    size = 2.45, linewidth = 0.3,
    color = ink, fill = "#E8F3EF"
  ) +
  annotate(
    "label", x = 5.65, y = 1,
    label = "Single-gene\ncomputational perturbation",
    size = 2.45, linewidth = 0.3,
    color = ink, fill = "#F4EEE0"
  ) +
  annotate(
    "label", x = 8.0, y = 1,
    label = "BH-adjusted downstream\nresponse-gene profiles",
    size = 2.45, linewidth = 0.3,
    color = ink, fill = "#F0EAF5"
  ) +
  annotate(
    "segment", x = 1.85, xend = 2.55, y = 1, yend = 1,
    arrow = arrow(length = unit(1.8, "mm")),
    color = mid_grey, linewidth = 0.45
  ) +
  annotate(
    "segment", x = 4.22, xend = 4.65, y = 1, yend = 1,
    arrow = arrow(length = unit(1.8, "mm")),
    color = mid_grey, linewidth = 0.45
  ) +
  annotate(
    "segment", x = 6.68, xend = 7.02, y = 1, yend = 1,
    arrow = arrow(length = unit(1.8, "mm")),
    color = mid_grey, linewidth = 0.45
  ) +
  coord_cartesian(xlim = c(0.1, 9), ylim = c(0.6, 1.4), clip = "off") +
  labs(title = "Analysis design") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(
      face = "bold", size = 8.6, color = ink,
      margin = margin(b = 1)
    ),
    plot.margin = margin(3, 5, 1, 5)
  )

panel_b <- ggplot(
  heatmap_dat,
  aes(x = candidate, y = gene, fill = pmax(pmin(Z, 4), -4))
) +
  geom_tile(color = "white", linewidth = 0.18) +
  geom_point(
    data = heatmap_dat[heatmap_dat$significant, , drop = FALSE],
    shape = 21, size = 0.95, stroke = 0.18,
    fill = ink, color = "white"
  ) +
  scale_fill_gradient2(
    low = "#356E9E", mid = "white", high = "#B4473F",
    midpoint = 0, limits = c(-4, 4),
    oob = squish, name = "Z score"
  ) +
  labs(
    title = "Downstream response-gene landscape",
    subtitle = paste0(
      "Display set: 12 strongest peak responses + 12 most variable responses;\n",
      "dots: BH-adjusted P < 0.05"
    ),
    x = NULL, y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(
      angle = 42, hjust = 1, vjust = 1, size = 6.1
    ),
    axis.text.y = element_text(size = 5.85),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(12, "mm")
  )

panel_c <- ggplot(
  rank_dat,
  aes(
    x = response_rank, y = neg_log10_q,
    group = candidate
  )
) +
  geom_hline(
    yintercept = -log10(0.05), linetype = "dashed",
    color = red, linewidth = 0.45
  ) +
  geom_line(
    data = rank_dat[rank_dat$candidate != "GSK3B", , drop = FALSE],
    color = "#98A4AE", linewidth = 0.42, alpha = 0.68
  ) +
  geom_line(
    data = rank_dat[rank_dat$candidate == "GSK3B", , drop = FALSE],
    color = target_colors[["GSK3B"]], linewidth = 1.0, alpha = 0.96
  ) +
  annotate(
    "text", x = 64, y = 1.52, label = "GSK3B",
    color = target_colors[["GSK3B"]],
    size = 2.25, hjust = 0
  ) +
  scale_x_continuous(
    breaks = c(1, 20, 40, 60, 80, 100, 120),
    limits = c(1, 120)
  ) +
  scale_y_continuous(
    limits = c(0, 12),
    breaks = c(0, -log10(0.05), 4, 8, 12),
    labels = c("0", "1.30", "4", "8", "\u226512")
  ) +
  labs(
    title = "Ranked response profiles",
    subtitle = "Dashed: BH-adjusted P = 0.05; gold: GSK3B",
    x = "Response-gene rank",
    y = expression(-log[10]("BH-adjusted P"))
  ) +
  theme_pub +
  theme(
    panel.grid.major.y = element_line(
      color = grid_col, linewidth = 0.28
    )
  )

panel_d <- ggplot(candidate_endpoints, aes(y = candidate)) +
  geom_segment(
    aes(
      x = 0,
      xend = log1p(candidate_endpoint),
      yend = candidate
    ),
    color = "#B9C1C8", linewidth = 0.75
  ) +
  geom_point(
    aes(
      x = log1p(candidate_endpoint),
      color = as.character(candidate)
    ),
    size = 2.6
  ) +
  geom_text(
    aes(
      x = log1p(candidate_endpoint),
      label = candidate_endpoint,
      color = as.character(candidate)
    ),
    hjust = -0.65, size = 2.15, show.legend = FALSE
  ) +
  scale_color_manual(values = target_colors, guide = "none") +
  scale_x_continuous(
    breaks = log1p(c(0, 5, 10, 20, 50, 100)),
    labels = c("0", "5", "10", "20", "50", "100"),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(
    title = "BH-significant downstream responses",
    subtitle = "All ten perturbations yielded significant response genes",
    x = "Response genes at BH-adjusted P < 0.05 (log1p axis)",
    y = NULL
  ) +
  theme_pub

panel_tag_theme <- theme(
  plot.tag = element_text(
    face = "bold", size = 8.8, color = ink, family = "sans"
  ),
  plot.tag.position = c(0.005, 0.995)
)
panel_a <- panel_a + labs(tag = "A") + panel_tag_theme
panel_b <- panel_b + labs(tag = "B") + panel_tag_theme
panel_c <- panel_c + labs(tag = "C") + panel_tag_theme
panel_d <- panel_d + labs(tag = "D") + panel_tag_theme

main_figure <- (
  panel_a /
    (panel_b | panel_c) /
    panel_d
) +
  plot_layout(heights = c(0.55, 2.9, 1.65)) +
  plot_annotation(
    title = "Computational single-gene perturbation reveals downstream network responses",
    subtitle = "Standard scTenifoldKnk analysis of ten prespecified candidates in the CD34\u207a HSPC network",
    theme = theme(
      plot.title = element_text(
        face = "bold", size = 12.8, color = ink,
        family = "sans"
      ),
      plot.subtitle = element_text(
        size = 7.7, color = muted, family = "sans"
      ),
      plot.tag = element_blank()
    )
  )

control_plot_a <- ggplot(
  control_plot,
  aes(x = log1p(response_count), y = candidate)
) +
  geom_boxplot(
    width = 0.55, outlier.shape = NA,
    fill = light_grey, color = mid_grey, linewidth = 0.35
  ) +
  geom_point(
    position = position_jitter(height = 0.12, width = 0),
    size = 0.72, alpha = 0.5, color = "#76828D"
  ) +
  geom_point(
    data = candidate_endpoints,
    aes(x = log1p(candidate_endpoint), y = candidate),
    inherit.aes = FALSE,
    shape = 23, size = 2.8, stroke = 0.5,
    fill = gold, color = ink
  ) +
  scale_x_continuous(
    breaks = log1p(c(0, 5, 10, 20, 50, 100)),
    labels = c("0", "5", "10", "20", "50", "100")
  ) +
  labs(
    title = "Candidate positions within assigned control groups",
    subtitle = "Diamonds: candidates; boxes and points: 20 controls per candidate",
    x = "Response genes (log1p axis)", y = NULL
  ) +
  theme_pub

supp_b <- ggplot(
  candidate_endpoints,
  aes(
    x = candidate_standardized_residual,
    y = candidate,
    color = as.character(candidate)
  )
) +
  geom_vline(
    xintercept = 0, color = "#9CA6AF", linewidth = 0.45
  ) +
  geom_segment(
    aes(x = 0, xend = candidate_standardized_residual, yend = candidate),
    linewidth = 0.68, alpha = 0.72
  ) +
  geom_point(
    aes(shape = individual_gate_pass),
    size = 2.5, stroke = 0.55
  ) +
  scale_color_manual(values = target_colors, guide = "none") +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(-0.7, 2.25)) +
  labs(
    title = "Standardized response residuals",
    subtitle = "TERT excluded from calibrated comparison",
    x = "Standardized residual", y = NULL
  ) +
  theme_pub

balance_before <- data.frame(
  candidate = original$candidate,
  imbalance = original$max_abs_standardized_matching_imbalance,
  assignment = "Initial assignment",
  stringsAsFactors = FALSE
)
balance_after <- data.frame(
  candidate = as.character(candidate_endpoints$candidate),
  imbalance = candidate_endpoints$max_abs_standardized_matching_imbalance,
  assignment = "Outcome-blind repartition",
  stringsAsFactors = FALSE
)
balance_plot <- bind_rows(balance_before, balance_after)
balance_plot$candidate <- factor(balance_plot$candidate, levels = rev(core10))
balance_plot$assignment <- factor(
  balance_plot$assignment,
  levels = c("Initial assignment", "Outcome-blind repartition")
)

supp_c <- ggplot(
  balance_plot,
  aes(x = imbalance, y = candidate, color = assignment)
) +
  geom_vline(
    xintercept = 1, linetype = "dashed",
    color = red, linewidth = 0.5
  ) +
  geom_line(
    aes(group = candidate),
    color = "#B9C1C8", linewidth = 0.42
  ) +
  geom_point(size = 2.05) +
  scale_color_manual(values = c(
    "Initial assignment" = "#9CA7B1",
    "Outcome-blind repartition" = teal
  )) +
  coord_cartesian(xlim = c(0, 2.9)) +
  labs(
    title = "Maximum baseline-covariate imbalance",
    subtitle = "Dashed line: absolute standardized difference of 1",
    x = "Maximum absolute standardized difference", y = NULL,
    color = NULL
  ) +
  theme_pub +
  theme(legend.position = "bottom")

supp_d <- ggplot(
  candidate_endpoints,
  aes(
    x = matched_control_local_percentile,
    y = candidate,
    color = ifelse(individual_gate_pass, "Calibration-eligible", "Not calibration-eligible")
  )
) +
  geom_vline(
    xintercept = 0.95, linetype = "dotted",
    color = red, linewidth = 0.5
  ) +
  geom_segment(
    aes(x = 0, xend = matched_control_local_percentile, yend = candidate),
    linewidth = 0.65, alpha = 0.72
  ) +
  geom_point(size = 2.4) +
  scale_color_manual(values = c(
    "Calibration-eligible" = blue,
    "Not calibration-eligible" = "#9CA7B1"
  )) +
  scale_x_continuous(
    limits = c(0, 1.02),
    breaks = c(0, 0.25, 0.5, 0.75, 0.95),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    title = "Local matched-control percentile",
    subtitle = "Dotted line marks the 95th percentile",
    x = "Local percentile", y = NULL, color = NULL
  ) +
  theme_pub +
  theme(legend.position = "bottom")

supp_e <- ggplot(
  balance_metric,
  aes(
    x = metric, y = candidate,
    fill = pmax(pmin(standardized_difference, 1.5), -1.5),
    color = abs(standardized_difference) > 1
  )
) +
  geom_tile(linewidth = 0.55) +
  geom_text(
    aes(label = sprintf("%.2f", standardized_difference)),
    size = 1.9, color = ink
  ) +
  scale_fill_gradient2(
    low = "#4079A8", mid = "white", high = "#BC584C",
    midpoint = 0, limits = c(-1.5, 1.5), oob = squish,
    name = "Standardized\ndifference"
  ) +
  scale_color_manual(
    values = c(`FALSE` = "white", `TRUE` = red),
    guide = "none"
  ) +
  labs(
    title = "Candidate-control balance for four matching variables",
    subtitle = "Red borders mark absolute standardized differences > 1",
    x = NULL, y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

null_plot <- pooled_null[
  pooled_null$eligible_group &
    is.finite(pooled_null$leave_one_out_standardized_residual),
  , drop = FALSE
]
calibrated_points <- sensitivity[
  sensitivity$individual_gate_pass, , drop = FALSE
]
supp_f <- ggplot(
  null_plot,
  aes(x = leave_one_out_standardized_residual)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 36, fill = light_grey, color = "white",
    linewidth = 0.25
  ) +
  geom_density(
    color = mid_grey, fill = "#B7C4CF",
    linewidth = 0.55, alpha = 0.28
  ) +
  geom_vline(
    xintercept = 0, color = "#9CA6AF", linewidth = 0.45
  ) +
  geom_rug(
    data = calibrated_points,
    aes(
      x = candidate_standardized_residual,
      color = candidate
    ),
    sides = "b", linewidth = 0.85, length = unit(3.5, "mm")
  ) +
  annotate(
    "text",
    x = sensitivity$candidate_standardized_residual[
      sensitivity$candidate == "GSK3B"
    ],
    y = Inf, label = "GSK3B",
    vjust = 1.45, hjust = -0.05,
    size = 2.1, color = target_colors[["GSK3B"]]
  ) +
  scale_color_manual(values = target_colors, guide = "none") +
  coord_cartesian(xlim = c(-1.5, 3.1)) +
  labs(
    title = "Pooled leave-one-out null distribution",
    subtitle = "Colored rug marks: nine calibrated candidate residuals",
    x = "Leave-one-out standardized residual",
    y = "Density"
  ) +
  theme_pub +
  theme(
    panel.grid.major.y = element_line(
      color = grid_col, linewidth = 0.28
    )
  )

supp_figure <- (
  (control_plot_a | supp_b) /
    (supp_c | supp_d) /
    (supp_e | supp_f)
) +
  plot_layout(heights = c(1.9, 1.75, 2.05)) +
  plot_annotation(
    title = "Matched-control sensitivity analysis of computational perturbation responses",
    subtitle = paste(
      "Ten candidates, 200 disjoint noncandidate controls,",
      "and four outcome-independent matching variables"
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(
        face = "bold", size = 12.8, color = ink,
        family = "sans"
      ),
      plot.subtitle = element_text(
        size = 7.7, color = muted, family = "sans"
      ),
      plot.tag = element_text(
        face = "bold", size = 10.8, color = ink,
        family = "sans"
      )
    )
  )

save_figure <- function(plot_object, stem, width_in, height_in) {
  svglite::svglite(
    paste0(stem, ".svg"),
    width = width_in, height = height_in, bg = "white"
  )
  print(plot_object)
  dev.off()

  grDevices::cairo_pdf(
    filename = paste0(stem, ".pdf"),
    width = width_in, height = height_in,
    family = "sans", onefile = TRUE
  )
  print(plot_object)
  dev.off()

  ragg::agg_png(
    paste0(stem, ".png"),
    width = width_in, height = height_in,
    units = "in", res = 500, background = "white"
  )
  print(plot_object)
  dev.off()

  ragg::agg_tiff(
    paste0(stem, ".tiff"),
    width = width_in, height = height_in,
    units = "in", res = 600, compression = "lzw",
    background = "white"
  )
  print(plot_object)
  dev.off()
}

save_figure(main_figure, main_stem, 7.205, 7.1)
save_figure(supp_figure, supp_stem, 7.205, 8.9)

cat("FIGURE8_STANDARD_AND_S11_SENSITIVITY_CANDIDATE_COMPLETE\n")
cat(paste(c(main_outputs, supp_outputs, source_paths), collapse = "\n"), "\n")
