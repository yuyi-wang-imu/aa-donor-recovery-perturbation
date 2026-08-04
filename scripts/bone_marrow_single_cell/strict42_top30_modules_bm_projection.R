# 93_strict42_Top30modules_BM_projection_20260716.R
#
# Purpose:
#   Recompute descriptive localization of the 11 strict42 WGCNA modules that
#   contain at least one formal Top30 candidate in the full GSE247531 bone-
#   marrow count matrix. This script writes a distinct set of 93-prefixed
#   outputs and stops if any output path already exists.
#
# Interpretation:
#   - Marker-defined coarse marrow compartments, not reference-mapped labels.
#   - Descriptive projection, not cell-type-specific WGCNA or causal evidence.
#   - Cells are not treated as independent subjects. A subject-by-timepoint
#     aggregation check is included for edge-release decisions.

suppressPackageStartupMessages({
  library(Matrix)
})

options(stringsAsFactors = FALSE)
set.seed(20260716)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript strict42_top30_modules_bm_projection.R <strict_wgcna_result_dir> <GSE247531_BMcounts.Rdata.gz>")
}
wgcna_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
bm_file <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
cmd_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_all, value = TRUE)
script_file <- if (length(file_arg) == 1L) normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE) else NA_character_
module_file <- file.path(wgcna_dir, "92_CD34_WGCNA_strict_donor_aware_20260715_module_gene_table.csv")
top30_file <- file.path(wgcna_dir, "92_CD34_WGCNA_strict_donor_aware_20260715_original_top30_strict_status.csv")
metadata_file <- file.path(wgcna_dir, "81_BM_sample_metadata_from_GEO.csv")

prefix <- file.path(wgcna_dir, "93_strict42_Top30modules_BM_projection_20260716")
output_files <- c(
  long = paste0(prefix, "_sensitivity_long.csv"),
  top1 = paste0(prefix, "_sensitivity_top1.csv"),
  agreement = paste0(prefix, "_agreement_summary.csv"),
  subject_timepoint = paste0(prefix, "_subject_timepoint_compartment_summary.csv"),
  hub_genes = paste0(prefix, "_module_hub_genes_used.csv"),
  celltype_summary = paste0(prefix, "_marker_celltype_summary.csv"),
  decisions = paste0(prefix, "_edge_decisions.csv"),
  formal_edges = paste0(prefix, "_formal_module_compartment_edges.csv"),
  source_hashes = paste0(prefix, "_source_hashes.tsv"),
  session_info = paste0(prefix, "_sessionInfo.txt"),
  run_log = paste0(prefix, "_run_log.txt"),
  output_hashes = paste0(prefix, "_output_hashes.tsv")
)

existing <- output_files[file.exists(output_files)]
if (length(existing) > 0) {
  stop("Refusing to overwrite existing 93-prefixed outputs: ", paste(existing, collapse = "; "))
}
for (f in c(script_file, module_file, top30_file, metadata_file, bm_file)) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}

log_con <- file(output_files[["run_log"]], open = "wt", encoding = "UTF-8")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
log_is_open <- TRUE
on.exit({
  if (isTRUE(log_is_open)) {
    try(sink(type = "message"), silent = TRUE)
    try(sink(), silent = TRUE)
    try(close(log_con), silent = TRUE)
  }
}, add = TRUE)

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")
}

write_tsv <- function(x, path) {
  write.table(
    x, path, sep = "\t", row.names = FALSE, quote = FALSE,
    fileEncoding = "UTF-8"
  )
}

cat("strict42 Top30-module BM projection started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R:", R.version.string, "\n")
cat("Seed: 20260716\n")
cat("Margin downgrade rule: absolute top1-top2 delta <0.01 OR relative delta <5%\n")
cat("Formal edge rule: >=8/9 top1 agreement, non-small baseline margin, and pooled top1 matching subject-timepoint median top1\n")

module_table <- read.csv(module_file, check.names = FALSE, fileEncoding = "UTF-8")
top30 <- read.csv(top30_file, check.names = FALSE, fileEncoding = "UTF-8")
bm_metadata <- read.csv(metadata_file, check.names = FALSE, fileEncoding = "UTF-8")

strict_modules <- sort(unique(top30$own_module_if_top5000[
  !is.na(top30$own_module_if_top5000) &
    nzchar(top30$own_module_if_top5000) &
    top30$own_module_if_top5000 != "NA"
]))
expected_modules <- sort(c(
  "blue", "darkred", "greenyellow", "lightcyan", "lightyellow", "magenta",
  "midnightblue", "purple", "royalblue", "turquoise", "yellow"
))
if (!identical(strict_modules, expected_modules)) {
  stop(
    "Unexpected strict module set. Observed: ", paste(strict_modules, collapse = ", "),
    "; expected: ", paste(expected_modules, collapse = ", ")
  )
}
cat("Strict Top30 modules (n=", length(strict_modules), "): ", paste(strict_modules, collapse = ", "), "\n", sep = "")

top_n_values <- c(30L, 50L, 100L)
get_module_hubs <- function(module_color, top_n) {
  kcol <- paste0("kME", module_color)
  d <- module_table[module_table$module_color == module_color, , drop = FALSE]
  if (!kcol %in% colnames(d)) stop("Missing kME column: ", kcol)
  kval <- suppressWarnings(as.numeric(d[[kcol]]))
  keep <- is.finite(kval) & !is.na(d$gene) & nzchar(d$gene)
  d <- d[keep, , drop = FALSE]
  kval <- kval[keep]
  ord <- order(abs(kval), decreasing = TRUE)
  unique(head(d$gene[ord], top_n))
}

module_hubs <- setNames(vector("list", length(top_n_values)), as.character(top_n_values))
for (top_n in top_n_values) {
  module_hubs[[as.character(top_n)]] <- setNames(
    lapply(strict_modules, get_module_hubs, top_n = top_n),
    strict_modules
  )
}

hub_records <- do.call(rbind, lapply(top_n_values, function(top_n) {
  do.call(rbind, lapply(strict_modules, function(module_color) {
    genes <- module_hubs[[as.character(top_n)]][[module_color]]
    data.frame(
      module_color = module_color,
      top_n_hubs = top_n,
      rank_within_selected_hubs = seq_along(genes),
      gene = genes,
      stringsAsFactors = FALSE
    )
  }))
}))

marker_sets <- list(
  HSPC = c("CD34", "KIT", "PROM1", "AVP", "GATA2", "MECOM", "MEIS1", "HLF", "SPINK2"),
  Erythroid = c("HBB", "HBA1", "HBA2", "ALAS2", "GYPA", "KLF1", "GATA1", "AHSP"),
  Megakaryocyte = c("PPBP", "PF4", "ITGA2B", "GP9", "MPL", "VWF"),
  Myeloid = c("LYZ", "LST1", "S100A8", "S100A9", "FCGR3A", "CD14", "MS4A7", "C1QA", "C1QB", "MPO", "ELANE", "AZU1", "CTSG", "FCER1G"),
  T_NK = c("CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", "NKG7", "GNLY", "PRF1", "GZMB", "KLRD1", "IL7R", "CCR7"),
  B_Plasma = c("MS4A1", "CD79A", "CD79B", "CD74", "BANK1", "MZB1", "JCHAIN", "IGHM", "IGKC", "SDC1", "XBP1"),
  Stromal_Endothelial = c("COL1A1", "COL1A2", "DCN", "LUM", "PECAM1", "VWF", "KDR", "CXCL12", "LEPR", "PDGFRA")
)
cell_type_order <- names(marker_sets)

selected_genes <- unique(c(
  unlist(marker_sets, use.names = FALSE),
  unlist(module_hubs[["100"]], use.names = FALSE)
))
cat("Selected scoring genes before BM matching:", length(selected_genes), "\n")

cat("Loading read-only BM sparse matrix:", bm_file, "\n")
load_env <- new.env(parent = emptyenv())
outer_con <- gzfile(bm_file, "rb")
inner_con <- gzcon(outer_con)
loaded_objects <- load(inner_con, envir = load_env)
close(inner_con)
if (!"counts" %in% loaded_objects) {
  stop("Expected object 'counts' not found. Loaded: ", paste(loaded_objects, collapse = ", "))
}
counts <- get("counts", envir = load_env)
rm(load_env)
gc()

cat("BM matrix class:", paste(class(counts), collapse = "/"), "\n")
cat("BM matrix dim:", paste(dim(counts), collapse = " x "), "\n")
cat("BM matrix nnzero:", Matrix::nnzero(counts), "\n")
cat("BM matrix object size GB:", round(as.numeric(object.size(counts)) / 1024^3, 3), "\n")

cell_ids <- colnames(counts)
count_suffix <- sub("^.*_([0-9]+)$", "_\\1", cell_ids)
count_suffix[!grepl("^_[0-9]+$", count_suffix)] <- NA_character_
metadata_index <- match(count_suffix, bm_metadata$count_suffix)
if (anyNA(metadata_index)) {
  bad <- unique(count_suffix[is.na(metadata_index)])
  stop("Unmapped BM count suffixes: ", paste(head(bad, 20), collapse = ", "))
}
cell_subject <- bm_metadata$subject[metadata_index]
cell_timepoint <- bm_metadata$timepoint[metadata_index]
cell_timepoint[is.na(cell_timepoint)] <- ""
subject_timepoint_key <- paste(cell_subject, cell_timepoint, sep = "||")

lib_size <- Matrix::colSums(counts)
lib_size[!is.finite(lib_size) | lib_size <= 0] <- 1
present_genes <- intersect(selected_genes, rownames(counts))
missing_genes <- setdiff(selected_genes, present_genes)
cat("Selected genes present:", length(present_genes), "; missing:", length(missing_genes), "\n")

expr <- counts[present_genes, , drop = FALSE]
rm(counts)
gc()

scale_by_cell <- 10000 / lib_size
expr@x <- log1p(expr@x * rep.int(scale_by_cell, diff(expr@p)))
rm(scale_by_cell, lib_size)
gc()

score_mean <- function(genes, mat = expr) {
  genes_present <- intersect(genes, rownames(mat))
  if (length(genes_present) < 1) return(rep(NA_real_, ncol(mat)))
  as.numeric(Matrix::colMeans(mat[genes_present, , drop = FALSE]))
}

score_detection <- function(genes, mat = expr) {
  genes_present <- intersect(genes, rownames(mat))
  if (length(genes_present) < 1) return(rep(NA_real_, ncol(mat)))
  as.numeric(Matrix::colMeans(mat[genes_present, , drop = FALSE] > 0))
}

marker_scores <- sapply(marker_sets, score_mean)
marker_scores <- as.matrix(marker_scores)
colnames(marker_scores) <- names(marker_sets)
best_idx <- max.col(marker_scores, ties.method = "first")
best_score <- marker_scores[cbind(seq_len(nrow(marker_scores)), best_idx)]
second_score <- apply(marker_scores, 1, function(x) sort(x, decreasing = TRUE)[2])
score_margin <- best_score - second_score
coarse_cell_type <- factor(colnames(marker_scores)[best_idx], levels = cell_type_order)

label_masks <- list(
  all_cells = rep(TRUE, length(coarse_cell_type)),
  confident_best_ge_0.02_margin_ge_0.005 = best_score >= 0.02 & score_margin >= 0.005,
  strict_best_ge_0.05_margin_ge_0.02 = best_score >= 0.05 & score_margin >= 0.02
)

celltype_summary <- do.call(rbind, lapply(cell_type_order, function(ct) {
  idx <- coarse_cell_type == ct
  data.frame(
    coarse_cell_type = ct,
    n_cells = sum(idx, na.rm = TRUE),
    pct_cells = 100 * mean(idx, na.rm = TRUE),
    n_subject_timepoints = length(unique(subject_timepoint_key[idx])),
    n_confident = sum(idx & label_masks[["confident_best_ge_0.02_margin_ge_0.005"]], na.rm = TRUE),
    n_strict = sum(idx & label_masks[["strict_best_ge_0.05_margin_ge_0.02"]], na.rm = TRUE),
    median_best_marker_score = ifelse(any(idx), median(best_score[idx], na.rm = TRUE), NA_real_),
    median_marker_margin = ifelse(any(idx), median(score_margin[idx], na.rm = TRUE), NA_real_),
    stringsAsFactors = FALSE
  )
}))
write_csv(celltype_summary, output_files[["celltype_summary"]])

for (top_n in top_n_values) {
  for (module_color in strict_modules) {
    genes <- module_hubs[[as.character(top_n)]][[module_color]]
    present <- intersect(genes, rownames(expr))
    idx <- hub_records$module_color == module_color & hub_records$top_n_hubs == top_n
    hub_records$present_in_BMcounts[idx] <- hub_records$gene[idx] %in% present
  }
}
write_csv(hub_records, output_files[["hub_genes"]])

summarize_scores <- function(
    score_vec, labels, keep, module_color, variant_id, sensitivity_family,
    top_n_hubs, scoring_method, label_strategy) {
  do.call(rbind, lapply(cell_type_order, function(ct) {
    idx <- keep & labels == ct & is.finite(score_vec)
    data.frame(
      module_color = module_color,
      variant_id = variant_id,
      sensitivity_family = sensitivity_family,
      top_n_hubs = top_n_hubs,
      scoring_method = scoring_method,
      label_strategy = label_strategy,
      coarse_cell_type = ct,
      n_cells = sum(idx, na.rm = TRUE),
      mean_score = ifelse(any(idx, na.rm = TRUE), mean(score_vec[idx], na.rm = TRUE), NA_real_),
      median_score = ifelse(any(idx, na.rm = TRUE), median(score_vec[idx], na.rm = TRUE), NA_real_),
      pct_cells_score_gt0 = ifelse(any(idx, na.rm = TRUE), 100 * mean(score_vec[idx] > 0, na.rm = TRUE), NA_real_),
      stringsAsFactors = FALSE
    )
  }))
}

mean_scores <- setNames(vector("list", length(top_n_values)), as.character(top_n_values))
for (top_n in top_n_values) {
  mean_scores[[as.character(top_n)]] <- setNames(
    lapply(strict_modules, function(module_color) {
      score_mean(module_hubs[[as.character(top_n)]][[module_color]])
    }),
    strict_modules
  )
}
detection_scores_50 <- setNames(
  lapply(strict_modules, function(module_color) {
    score_detection(module_hubs[["50"]][[module_color]])
  }),
  strict_modules
)

records <- list()
add_record <- function(x) {
  records[[length(records) + 1L]] <<- x
}

cat("Computing three hub-count variants...\n")
for (top_n in top_n_values) {
  for (module_color in strict_modules) {
    add_record(summarize_scores(
      mean_scores[[as.character(top_n)]][[module_color]],
      coarse_cell_type, label_masks[["all_cells"]], module_color,
      paste0("hub_top", top_n, "_mean_all"), "hub_count", top_n,
      "mean_log1pCP10K", "all_cells"
    ))
  }
}

cat("Computing two scoring-method variants...\n")
for (module_color in strict_modules) {
  add_record(summarize_scores(
    mean_scores[["50"]][[module_color]], coarse_cell_type,
    label_masks[["all_cells"]], module_color,
    "scoring_top50_mean_all", "scoring_method", 50,
    "mean_log1pCP10K", "all_cells"
  ))
  add_record(summarize_scores(
    detection_scores_50[[module_color]], coarse_cell_type,
    label_masks[["all_cells"]], module_color,
    "scoring_top50_detection_all", "scoring_method", 50,
    "detection_fraction", "all_cells"
  ))
}

cat("Computing three label-threshold variants...\n")
label_variant_ids <- c(
  all_cells = "label_top50_mean_all",
  confident_best_ge_0.02_margin_ge_0.005 = "label_top50_mean_confident",
  strict_best_ge_0.05_margin_ge_0.02 = "label_top50_mean_strict"
)
for (label_strategy in names(label_masks)) {
  for (module_color in strict_modules) {
    add_record(summarize_scores(
      mean_scores[["50"]][[module_color]], coarse_cell_type,
      label_masks[[label_strategy]], module_color,
      label_variant_ids[[label_strategy]], "label_threshold", 50,
      "mean_log1pCP10K", label_strategy
    ))
  }
}

cat("Computing sampled-rank variant...\n")
sample_per_type <- 3000L
sample_idx <- unlist(lapply(cell_type_order, function(ct) {
  idx <- which(coarse_cell_type == ct)
  if (length(idx) <= sample_per_type) idx else sample(idx, sample_per_type)
}), use.names = FALSE)
sample_idx <- sort(unique(sample_idx))
sample_labels <- coarse_cell_type[sample_idx]
sample_mat <- as.matrix(expr[, sample_idx, drop = FALSE])
rank_mat <- apply(sample_mat, 2, function(v) rank(v, ties.method = "average") / length(v))
rownames(rank_mat) <- rownames(sample_mat)
rm(sample_mat)
gc()

rank_score <- function(genes) {
  genes_present <- intersect(genes, rownames(rank_mat))
  if (length(genes_present) < 1) return(rep(NA_real_, ncol(rank_mat)))
  colMeans(rank_mat[genes_present, , drop = FALSE], na.rm = TRUE)
}
for (module_color in strict_modules) {
  add_record(summarize_scores(
    rank_score(module_hubs[["50"]][[module_color]]), sample_labels,
    rep(TRUE, length(sample_labels)), module_color,
    "rank_top50_sampled_all", "sampled_rank", 50,
    "sampled_rank_percentile", "sampled_all_cells"
  ))
}

sensitivity_long <- do.call(rbind, records)
sensitivity_long <- sensitivity_long[order(
  sensitivity_long$module_color,
  sensitivity_long$variant_id,
  match(sensitivity_long$coarse_cell_type, cell_type_order)
), ]
if (length(unique(sensitivity_long$variant_id)) != 9L) {
  stop("Expected exactly 9 sensitivity variants; observed ", length(unique(sensitivity_long$variant_id)))
}
write_csv(sensitivity_long, output_files[["long"]])

top1_records <- list()
for (module_color in strict_modules) {
  for (variant_id in unique(sensitivity_long$variant_id)) {
    d <- sensitivity_long[
      sensitivity_long$module_color == module_color &
        sensitivity_long$variant_id == variant_id,
      , drop = FALSE
    ]
    d <- d[order(d$mean_score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    top1_records[[length(top1_records) + 1L]] <- data.frame(
      module_color = module_color,
      variant_id = variant_id,
      sensitivity_family = d$sensitivity_family[1],
      top_n_hubs = d$top_n_hubs[1],
      scoring_method = d$scoring_method[1],
      label_strategy = d$label_strategy[1],
      top1_cell_type = d$coarse_cell_type[1],
      top1_mean_score = d$mean_score[1],
      top2_cell_type = d$coarse_cell_type[2],
      top2_mean_score = d$mean_score[2],
      top1_top2_delta = d$mean_score[1] - d$mean_score[2],
      top1_top2_relative_delta = (d$mean_score[1] - d$mean_score[2]) / max(abs(d$mean_score[1]), .Machine$double.eps),
      stringsAsFactors = FALSE
    )
  }
}
top1 <- do.call(rbind, top1_records)
baseline_variant <- "hub_top50_mean_all"
base <- top1[top1$variant_id == baseline_variant, c(
  "module_color", "top1_cell_type", "top1_mean_score", "top2_cell_type",
  "top2_mean_score", "top1_top2_delta", "top1_top2_relative_delta"
)]
names(base)[-1] <- paste0("baseline_", names(base)[-1])
top1 <- merge(top1, base, by = "module_color", all.x = TRUE, sort = FALSE)
top1$top1_agrees_with_baseline <- top1$top1_cell_type == top1$baseline_top1_cell_type
top1 <- top1[order(top1$module_color, top1$variant_id), ]
write_csv(top1, output_files[["top1"]])

cat("Computing subject-by-timepoint baseline localization check...\n")
subject_timepoint_records <- list()
subject_timepoint_top <- list()
group_key <- paste(subject_timepoint_key, as.character(coarse_cell_type), sep = "@@")
for (module_color in strict_modules) {
  x <- mean_scores[["50"]][[module_color]]
  rs <- rowsum(cbind(score_sum = x, n_cells = rep.int(1, length(x))), group_key, reorder = FALSE)
  agg <- data.frame(group_key = rownames(rs), rs, stringsAsFactors = FALSE)
  parts <- strsplit(agg$group_key, "@@", fixed = TRUE)
  agg$subject_timepoint <- vapply(parts, `[`, character(1), 1)
  agg$coarse_cell_type <- vapply(parts, `[`, character(1), 2)
  agg$mean_score <- agg$score_sum / agg$n_cells
  agg$module_color <- module_color
  subject_timepoint_records[[length(subject_timepoint_records) + 1L]] <- agg[, c(
    "module_color", "subject_timepoint", "coarse_cell_type", "n_cells", "mean_score"
  )]

  med <- do.call(rbind, lapply(cell_type_order, function(ct) {
    z <- agg$mean_score[agg$coarse_cell_type == ct & is.finite(agg$mean_score)]
    data.frame(
      coarse_cell_type = ct,
      n_subject_timepoints = length(z),
      median_subject_timepoint_mean_score = ifelse(length(z) > 0, median(z), NA_real_),
      stringsAsFactors = FALSE
    )
  }))
  med <- med[order(med$median_subject_timepoint_mean_score, decreasing = TRUE, na.last = TRUE), ]
  subject_timepoint_top[[length(subject_timepoint_top) + 1L]] <- data.frame(
    module_color = module_color,
    subject_timepoint_top1_cell_type = med$coarse_cell_type[1],
    subject_timepoint_top1_median_score = med$median_subject_timepoint_mean_score[1],
    subject_timepoint_top2_cell_type = med$coarse_cell_type[2],
    subject_timepoint_top2_median_score = med$median_subject_timepoint_mean_score[2],
    subject_timepoint_top1_top2_delta = med$median_subject_timepoint_mean_score[1] - med$median_subject_timepoint_mean_score[2],
    stringsAsFactors = FALSE
  )
}
subject_timepoint_long <- do.call(rbind, subject_timepoint_records)
write_csv(subject_timepoint_long, output_files[["subject_timepoint"]])
subject_top <- do.call(rbind, subject_timepoint_top)

agreement <- do.call(rbind, lapply(strict_modules, function(module_color) {
  d <- top1[top1$module_color == module_color, , drop = FALSE]
  b <- d[d$variant_id == baseline_variant, , drop = FALSE][1, ]
  data.frame(
    module_color = module_color,
    baseline_top1_cell_type = b$baseline_top1_cell_type,
    baseline_top1_mean_score = b$baseline_top1_mean_score,
    baseline_top2_cell_type = b$baseline_top2_cell_type,
    baseline_top2_mean_score = b$baseline_top2_mean_score,
    baseline_top1_top2_delta = b$baseline_top1_top2_delta,
    baseline_top1_top2_relative_delta = b$baseline_top1_top2_relative_delta,
    n_variants_tested = nrow(d),
    n_variants_agree_top1 = sum(d$top1_agrees_with_baseline, na.rm = TRUE),
    pct_variants_agree_top1 = 100 * mean(d$top1_agrees_with_baseline, na.rm = TRUE),
    discordant_variants = paste(d$variant_id[!d$top1_agrees_with_baseline], collapse = "; "),
    stringsAsFactors = FALSE
  )
}))
agreement <- merge(agreement, subject_top, by = "module_color", all.x = TRUE, sort = FALSE)
agreement$agreement_ge_8_of_9 <- agreement$n_variants_agree_top1 >= 8
agreement$small_baseline_margin <-
  agreement$baseline_top1_top2_delta < 0.01 |
  agreement$baseline_top1_top2_relative_delta < 0.05
agreement$subject_timepoint_top1_matches_pooled <-
  agreement$subject_timepoint_top1_cell_type == agreement$baseline_top1_cell_type
agreement$edge_decision <- ifelse(
  !agreement$agreement_ge_8_of_9,
  "exclude_instability_below_8_of_9",
  ifelse(
    agreement$small_baseline_margin,
    "downgrade_small_top1_top2_margin",
    ifelse(
      !agreement$subject_timepoint_top1_matches_pooled,
      "downgrade_subject_timepoint_top1_discordance",
      "include_descriptive_context_edge"
    )
  )
)
agreement$is_formal_edge <- agreement$edge_decision == "include_descriptive_context_edge"
agreement <- agreement[order(agreement$module_color), ]
write_csv(agreement, output_files[["agreement"]])
write_csv(agreement, output_files[["decisions"]])

formal <- agreement[agreement$is_formal_edge, , drop = FALSE]
formal_edges <- data.frame(
  source = paste0("MOD_", formal$module_color),
  target = paste0("CELL_", formal$baseline_top1_cell_type),
  edge_type = "Module-compartment descriptive context",
  interaction = "marker-defined marrow localization",
  module_label = paste0("ME", formal$module_color),
  compartment_label = formal$baseline_top1_cell_type,
  n_variants_agree_top1 = formal$n_variants_agree_top1,
  n_variants_tested = formal$n_variants_tested,
  pct_variants_agree_top1 = formal$pct_variants_agree_top1,
  baseline_top1_mean_score = formal$baseline_top1_mean_score,
  baseline_top2_cell_type = formal$baseline_top2_cell_type,
  baseline_top2_mean_score = formal$baseline_top2_mean_score,
  baseline_top1_top2_delta = formal$baseline_top1_top2_delta,
  baseline_top1_top2_relative_delta = formal$baseline_top1_top2_relative_delta,
  subject_timepoint_top1_matches_pooled = formal$subject_timepoint_top1_matches_pooled,
  evidence_boundary = "Descriptive marker-defined compartment context; not cell-type-specific WGCNA, causal localization, or independent validation",
  stringsAsFactors = FALSE
)
write_csv(formal_edges, output_files[["formal_edges"]])

capture.output(sessionInfo(), file = output_files[["session_info"]])

cat("Hashing four read-only source inputs and the new script...\n")
source_paths <- c(script_file, module_file, top30_file, metadata_file, bm_file)
source_roles <- c(
  "analysis_script", "strict42_module_gene_table", "strict42_Top30_status",
  "BM_sample_metadata", "read_only_BM_count_matrix"
)
source_hash_values <- unname(tools::sha256sum(source_paths))
source_hashes <- data.frame(
  role = source_roles,
  file_name = basename(source_paths),
  bytes = file.info(source_paths)$size,
  sha256 = source_hash_values,
  stringsAsFactors = FALSE
)
write_tsv(source_hashes, output_files[["source_hashes"]])

cat("\nFormal edge decisions:\n")
print(agreement[, c(
  "module_color", "baseline_top1_cell_type", "n_variants_agree_top1",
  "baseline_top1_top2_delta", "baseline_top1_top2_relative_delta",
  "subject_timepoint_top1_cell_type", "edge_decision"
)], row.names = FALSE)
cat("Formal descriptive context edges released:", nrow(formal_edges), "\n")
cat("strict42 Top30-module BM projection finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

sink(type = "message")
sink()
close(log_con)
log_is_open <- FALSE

hash_targets <- unname(output_files[setdiff(names(output_files), "output_hashes")])
output_hash_table <- data.frame(
  file_name = basename(hash_targets),
  bytes = file.info(hash_targets)$size,
  sha256 = unname(tools::sha256sum(hash_targets)),
  stringsAsFactors = FALSE
)
write_tsv(output_hash_table, output_files[["output_hashes"]])
