#!/usr/bin/env Rscript

# Publication figure for the 12 two-herb antecedent association rules.
# All plotted quantities are read from the source-backed SPSS extraction JSON.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(jsonlite)
  library(scales)
  library(svglite)
  library(ragg)
  library(png)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript plot_verified_spss_association_rules.R <story_root>")
}
story_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
source_json <- file.path(story_dir, "05_tables", "association_rules_SPSS_verified_source_20260717_v1.json")
figure_dir <- file.path(story_dir, "03_main_figure_candidates")
figure_version <- Sys.getenv("FIGURE_VERSION", unset = "v1")
if (!grepl("^v[0-9]+$", figure_version)) stop("FIGURE_VERSION must look like v1, v2, ...")
output_base <- file.path(figure_dir, paste0("Figure_1C_association_rules_SPSS_source_R_20260717_", figure_version))
qa_path <- file.path(story_dir, "09_quality_control", paste0("Figure_1C_association_rules_SPSS_source_R_QA_20260717_", figure_version, ".json"))

required_outputs <- c(
  paste0(output_base, ".png"),
  paste0(output_base, ".tiff"),
  paste0(output_base, ".pdf"),
  paste0(output_base, ".svg"),
  qa_path
)
existing_outputs <- required_outputs[file.exists(required_outputs)]
if (length(existing_outputs) > 0) {
  stop(
    "Refusing to overwrite existing output(s): ",
    paste(existing_outputs, collapse = "; ")
  )
}

src <- fromJSON(source_json, simplifyDataFrame = TRUE)
rules <- as.data.frame(src$rules, stringsAsFactors = FALSE)

expected_abbrev <- c(
  "EH + CS \u2192 LLF",
  "CCC + AMR \u2192 ASR",
  "ACC + RRP \u2192 ASR",
  "AH + RRP \u2192 ASR",
  "ACC + CS \u2192 ASR",
  "MOR + RRP \u2192 ASR",
  "MOR + EF \u2192 ASR",
  "CR + POR \u2192 AMR",
  "EH + ACC \u2192 ASR",
  "EF + ACC \u2192 AR",
  "PRA + RRP \u2192 ASR",
  "PR + RRP \u2192 ASR"
)

expected_antecedent_n <- c(45, 42, 64, 40, 60, 49, 39, 48, 47, 47, 51, 42)
expected_joint_n <- c(44, 39, 58, 36, 54, 44, 35, 43, 42, 42, 45, 37)
expected_coverage <- c(
  11.5384615384615, 10.7692307692307, 16.4102564102564,
  10.2564102564102, 15.3846153846153, 12.5641025641025,
  10.0, 12.3076923076923, 12.051282051282,
  12.051282051282, 13.076923076923, 10.7692307692307
)
expected_support <- c(
  11.2820512820512, 10.0, 14.8717948717948,
  9.23076923076923, 13.8461538461538, 11.2820512820512,
  8.97435897435897, 11.025641025641, 10.7692307692307,
  10.7692307692307, 11.5384615384615, 9.48717948717948
)
expected_confidence <- c(
  97.7777777777777, 92.8571428571428, 90.625, 90.0, 90.0,
  89.7959183673469, 89.7435897435897, 89.5833333333333,
  89.3617021276595, 89.3617021276595, 88.235294117647,
  88.095238095238
)
expected_lift <- c(
  3.02645502645502, 1.51524208009563, 1.47881799163179,
  1.46861924686192, 1.46861924686192, 1.46528904448808,
  1.46443514644351, 2.3767006802721, 1.45820350752247,
  1.52188051658459, 1.4398227910411, 1.43753735803945
)

rule_abbrev <- c(
  "Ecliptae Herba" = "EH",
  "Cuscutae Semen" = "CS",
  "Ligustri Lucidi Fructus" = "LLF",
  "Cervi Cornus Colla" = "CCC",
  "Atractylodis Macrocephalae Rhizoma" = "AMR",
  "Angelicae Sinensis Radix" = "ASR",
  "Asini Corii Colla" = "ACC",
  "Rehmanniae Radix Praeparata" = "RRP",
  "Agrimoniae Herba" = "AH",
  "Morindae Officinalis Radix" = "MOR",
  "Epimedii Folium" = "EF",
  "Codonopsis Radix" = "CR",
  "Poria" = "POR",
  "Astragali Radix" = "AR",
  "Paeoniae Radix Alba" = "PRA",
  "Polygonati Rhizoma" = "PR"
)

get_abbrev <- function(x) {
  out <- unname(rule_abbrev[x])
  if (anyNA(out)) stop("An herb name lacks a frozen abbreviation: ", paste(x[is.na(out)], collapse = ", "))
  out
}

rules <- rules %>%
  mutate(
    rule_abbrev = paste0(
      get_abbrev(antecedent_1_en), " + ",
      get_abbrev(antecedent_2_en), " \u2192 ",
      get_abbrev(consequent_en)
    ),
    rule_label = paste0("R", rank, "   ", rule_abbrev),
    rule_label = factor(rule_label, levels = rev(rule_label)),
    is_core = rank == 1,
    confidence_label = sprintf("%.3f%%", confidence_pct)
  )

assert_close <- function(actual, expected, label, tol = 1e-10) {
  if (length(actual) != length(expected) || any(abs(actual - expected) > tol)) {
    stop(label, " does not match the frozen SPSS extraction.")
  }
}

if (nrow(rules) != 12L) stop("Expected 12 rules; observed ", nrow(rules), ".")
if (!identical(as.integer(rules$rank), 1:12)) stop("Rule ranks are not the frozen 1-12 order.")
if (!identical(rules$rule_abbrev, expected_abbrev)) stop("Rule identities/order differ from the frozen SPSS extraction.")
if (!all(rules$antecedent_1_en != "" & rules$antecedent_2_en != "")) stop("Every retained rule must have exactly two antecedent herbs.")
if (!all(rules$antecedent_coverage_pct >= 10)) stop("At least one rule falls below the 10% antecedent coverage threshold.")
if (!all(rules$confidence_pct >= 88)) stop("At least one rule falls below the 88% confidence threshold.")
assert_close(rules$antecedent_n, expected_antecedent_n, "Antecedent counts")
assert_close(rules$joint_n, expected_joint_n, "Joint counts")
assert_close(rules$antecedent_coverage_pct, expected_coverage, "Antecedent coverage")
assert_close(rules$joint_support_pct, expected_support, "Joint support")
assert_close(rules$confidence_pct, expected_confidence, "Confidence")
assert_close(rules$lift, expected_lift, "Lift")

palette <- c(
  ink = "#24313D",
  neutral = "#B8C2CC",
  grid = "#E7ECF0",
  blue = "#3E6F9E",
  teal = "#49A7A0",
  gold = "#E6A11A",
  gold_dark = "#9A5D00"
)

base_family <- "Arial"

p <- ggplot(rules, aes(x = confidence_pct, y = rule_label)) +
  geom_vline(
    xintercept = 88,
    linetype = "22",
    linewidth = 0.45,
    colour = "#7D8994"
  ) +
  geom_segment(
    aes(x = 88, xend = confidence_pct, yend = rule_label),
    linewidth = 0.75,
    lineend = "round",
    colour = palette[["neutral"]]
  ) +
  geom_point(
    aes(size = joint_support_pct, fill = lift),
    shape = 21,
    colour = "white",
    stroke = 0.7
  ) +
  geom_point(
    data = filter(rules, is_core),
    aes(size = joint_support_pct, fill = lift),
    shape = 21,
    colour = palette[["gold_dark"]],
    stroke = 1.25,
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = confidence_pct + 0.28, label = confidence_label),
    hjust = 0,
    size = 2.65,
    family = base_family,
    colour = palette[["ink"]]
  ) +
  scale_x_continuous(
    name = "Confidence (%)",
    breaks = c(88, 90, 92, 94, 96, 98, 100),
    limits = c(87.55, 100.55),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(name = NULL, expand = expansion(add = c(0.6, 0.65))) +
  scale_size_area(
    name = "Joint support (%)",
    max_size = 7.2,
    breaks = c(9, 11, 13, 15),
    limits = c(8.8, 15.1)
  ) +
  scale_fill_gradientn(
    name = "Lift",
    colours = c(palette[["blue"]], palette[["teal"]], palette[["gold"]]),
    values = rescale(c(min(rules$lift), 2.0, max(rules$lift))),
    limits = range(rules$lift),
    breaks = c(1.5, 2.0, 2.5, 3.0),
    oob = squish
  ) +
  guides(
    size = guide_legend(
      order = 1,
      title.position = "top",
      nrow = 1,
      override.aes = list(fill = "#8FC7C1", colour = "white", stroke = 0.6)
    ),
    fill = guide_colourbar(
      order = 2,
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(33, "mm"),
      barheight = unit(3.2, "mm"),
      ticks.colour = palette[["ink"]],
      frame.colour = "#AAB4BE"
    )
  ) +
  labs(
    title = "High-confidence two-herb association rules identified by Apriori",
    subtitle = paste0(
      "n = 390 prescription records; point area = joint support; color = lift; dashed line = 88% confidence threshold\n",
      "Core rule R1: n = 45/44; antecedent coverage = 11.538%; joint support = 11.282%; confidence = 97.778%; lift = 3.026"
    ),
    caption = "Rules are shown in the original SPSS confidence-ranked order. Herb-name abbreviations are defined in the figure legend."
  ) +
  theme_minimal(base_size = 8.2, base_family = base_family) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = palette[["grid"]], linewidth = 0.35),
    axis.title.x = element_text(size = 8.2, margin = margin(t = 6)),
    axis.text.x = element_text(size = 7.2, colour = palette[["ink"]]),
    axis.text.y = element_text(size = 7.1, colour = palette[["ink"]], hjust = 1),
    axis.ticks.x = element_line(linewidth = 0.35, colour = palette[["ink"]]),
    axis.ticks.length = unit(1.5, "mm"),
    plot.title = element_text(size = 10.5, face = "bold", colour = palette[["ink"]], margin = margin(b = 4)),
    plot.subtitle = element_text(size = 7.3, colour = "#52606C", lineheight = 1.18, margin = margin(b = 8)),
    plot.caption = element_text(size = 6.5, colour = "#65727E", hjust = 0, margin = margin(t = 7)),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.title = element_text(size = 7.1, face = "bold", colour = palette[["ink"]]),
    legend.text = element_text(size = 6.7, colour = palette[["ink"]]),
    legend.key = element_blank(),
    legend.spacing.x = unit(2.5, "mm"),
    legend.margin = margin(t = 2, r = 0, b = 0, l = 0),
    plot.margin = margin(t = 9, r = 10, b = 7, l = 7, unit = "mm")
  )

width_mm <- 183
height_mm <- 130
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
dpi <- 600

ragg::agg_png(
  paste0(output_base, ".png"),
  width = width_in,
  height = height_in,
  units = "in",
  res = dpi,
  background = "white",
  scaling = 1
)
print(p)
dev.off()

ragg::agg_tiff(
  paste0(output_base, ".tiff"),
  width = width_in,
  height = height_in,
  units = "in",
  res = dpi,
  compression = "lzw",
  background = "white",
  scaling = 1
)
print(p)
dev.off()

svglite::svglite(
  paste0(output_base, ".svg"),
  width = width_in,
  height = height_in,
  bg = "white"
)
print(p)
dev.off()

grDevices::cairo_pdf(
  paste0(output_base, ".pdf"),
  width = width_in,
  height = height_in,
  family = base_family,
  bg = "white",
  onefile = TRUE
)
print(p)
dev.off()

expected_pixels <- c(width = floor(width_in * dpi), height = floor(height_in * dpi))
png_native <- png::readPNG(paste0(output_base, ".png"), native = TRUE, info = TRUE)
png_dimensions <- dim(png_native)
actual_pixels <- c(width = png_dimensions[2], height = png_dimensions[1])

qa <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  backend = R.version.string,
  source_file = basename(source_json),
  source_sha256 = src$source_sha256,
  source_output_rule_count = src$source_output_rule_count,
  selected_rule_count = nrow(rules),
  rule_order_exact = identical(as.integer(rules$rank), 1:12),
  rule_identity_exact = identical(rules$rule_abbrev, expected_abbrev),
  all_exactly_two_antecedents = all(rules$antecedent_1_en != "" & rules$antecedent_2_en != ""),
  all_coverage_ge_10 = all(rules$antecedent_coverage_pct >= 10),
  all_confidence_ge_88 = all(rules$confidence_pct >= 88),
  core_rule = list(
    rank = rules$rank[1],
    rule = rules$rule_abbrev[1],
    antecedent_n = rules$antecedent_n[1],
    joint_n = rules$joint_n[1],
    antecedent_coverage_pct = rules$antecedent_coverage_pct[1],
    joint_support_pct = rules$joint_support_pct[1],
    confidence_pct = rules$confidence_pct[1],
    lift = rules$lift[1]
  ),
  visual_encodings = list(
    x = "confidence_pct",
    point_area = "joint_support_pct",
    point_fill = "lift",
    vertical_dashed_line = "88% confidence threshold",
    core_rule = "gold-dark outline plus a metrics line above the plotting region"
  ),
  export = list(
    width_mm = width_mm,
    height_mm = height_mm,
    dpi = dpi,
    expected_png_pixels = as.list(expected_pixels),
    actual_png_pixels = as.list(actual_pixels),
    png_dimensions_match = identical(unname(as.integer(actual_pixels)), unname(as.integer(expected_pixels))),
    files = lapply(required_outputs[1:4], function(x) {
      list(file = basename(x), exists = file.exists(x), bytes = unname(file.info(x)$size))
    })
  )
)

write_json(qa, qa_path, pretty = TRUE, auto_unbox = TRUE, digits = 15)

cat("Generated:\n", paste(required_outputs, collapse = "\n"), "\n", sep = "")
