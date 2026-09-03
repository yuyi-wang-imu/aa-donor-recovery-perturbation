invisible(Sys.setlocale("LC_ALL", "Chinese"))

suppressPackageStartupMessages({
  library(WGCNA)
})

options(stringsAsFactors = FALSE)
set.seed(20260902)
disableWGCNAThreads()
allowWGCNAThreads(nThreads = 4)

root <- "E:/AA\u65b0\u6295\u7a3f\u601d\u8def"
qc_dir <- file.path(root, "BMC_Pharmacology_and_Toxicology_\u6295\u7a3f_20260830", "90_\u5185\u90e8QC_\u52ff\u63d0\u4ea4")
hspc_dir <- file.path(
  root,
  "02_\u9879\u76ee\u603b\u5f52\u6863_\u62cd\u677f\u7248",
  "04_\u524d\u7f6e\u65b9\u6cd5\u4e0e\u7b97\u6cd5",
  "\u65b9\u6848\u524d\u7f6e\u65b9\u6cd5",
  "02_GEO_WGCNA_\u539f\u59cb\u6570\u636e\u4e0e\u8bbe\u8ba1",
  "04_GSE247531_CD34_HSPC\u5206\u6790\u7ed3\u679c"
)

prefix <- "simulated_editor_baseline_only_wgcna_20260902_v1"
meta_path <- file.path(hspc_dir, "92_CD34_WGCNA_strict_donor_aware_20260715_collapsed_sample_metadata.tsv")
counts_path <- file.path(hspc_dir, "92_CD34_WGCNA_strict_donor_aware_20260715_collapsed_counts_by_subject_timepoint.csv")
original_module_path <- file.path(hspc_dir, "92_CD34_WGCNA_strict_donor_aware_20260715_module_gene_table.csv")

outputs <- c(
  file.path(qc_dir, paste0(prefix, "_soft_threshold_fit.csv")),
  file.path(qc_dir, paste0(prefix, "_module_gene_table.csv")),
  file.path(qc_dir, paste0(prefix, "_module_sizes.csv")),
  file.path(qc_dir, paste0(prefix, "_original42_overlap_best_pairs.csv")),
  file.path(qc_dir, paste0(prefix, "_candidate10_mapping.csv")),
  file.path(qc_dir, paste0(prefix, "_summary.tsv")),
  file.path(qc_dir, paste0(prefix, "_report.md")),
  file.path(qc_dir, paste0(prefix, "_runtime.txt"))
)
existing <- outputs[file.exists(outputs)]
if (length(existing) > 0) {
  stop("Refusing to overwrite existing outputs: ", paste(existing, collapse = "; "))
}

started <- Sys.time()
meta <- read.delim(meta_path, check.names = FALSE)
counts_df <- read.csv(counts_path, check.names = FALSE)
gene <- toupper(as.character(counts_df[[1]]))
counts <- as.matrix(counts_df[, -1, drop = FALSE])
storage.mode(counts) <- "numeric"
rownames(counts) <- gene
if (any(duplicated(rownames(counts)))) {
  counts <- rowsum(counts, group = rownames(counts), reorder = FALSE)
}

keep_meta <- meta$disease == "HD" | (meta$disease == "SAA" & meta$timepoint == "baseline")
baseline_meta <- meta[keep_meta, , drop = FALSE]
sample_ids <- as.character(baseline_meta$subject_timepoint_id)
if (!all(sample_ids %in% colnames(counts))) {
  stop("Collapsed count matrix does not contain all baseline/reference profiles")
}
counts <- counts[, sample_ids, drop = FALSE]

if (nrow(baseline_meta) != length(unique(baseline_meta$subject))) {
  stop("Baseline-only subset is not one profile per participant")
}
if (sum(baseline_meta$disease == "HD") != 4 || sum(baseline_meta$disease == "SAA") != 19) {
  stop("Unexpected baseline/reference group sizes")
}

lib_size <- colSums(counts)
log1p_cp10k <- log1p(sweep(counts, 2, lib_size / 10000, "/"))
log1p_cp10k[!is.finite(log1p_cp10k)] <- 0

detected_samples <- rowSums(counts > 0)
total_counts <- rowSums(counts)
is_mito <- grepl("^MT-", rownames(counts))
keep_expr <- detected_samples >= 8 & total_counts >= 50 & !is_mito
mad_value <- apply(log1p_cp10k[keep_expr, , drop = FALSE], 1, mad, na.rm = TRUE)
mad_value[!is.finite(mad_value)] <- 0
wgcna_genes <- names(sort(mad_value, decreasing = TRUE))[seq_len(min(5000, length(mad_value)))]

datExpr <- as.data.frame(t(log1p_cp10k[wgcna_genes, sample_ids, drop = FALSE]), check.names = FALSE)
rownames(datExpr) <- sample_ids
gsg <- goodSamplesGenes(datExpr, verbose = 1)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}

powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = list(maxPOutliers = 0.1, use = "pairwise.complete.obs"),
  verbose = 2
)
fit_indices <- sft$fitIndices
fit_indices$fixed_power_used <- 18
fit_indices$reestimated_powerEstimate <- sft$powerEstimate
write.csv(fit_indices, outputs[1], row.names = FALSE, quote = TRUE)

net <- blockwiseModules(
  datExpr,
  power = 18,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  maxPOutliers = 0.1,
  maxBlockSize = ncol(datExpr),
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 2
)
module_colors <- net$colors
names(module_colors) <- colnames(datExpr)
module_table <- data.frame(gene = names(module_colors), baseline_module = unname(module_colors))
write.csv(module_table, outputs[2], row.names = FALSE, quote = TRUE)

module_sizes <- as.data.frame(table(module_table$baseline_module), stringsAsFactors = FALSE)
names(module_sizes) <- c("baseline_module", "n_genes")
module_sizes <- module_sizes[order(-module_sizes$n_genes, module_sizes$baseline_module), ]
write.csv(module_sizes, outputs[3], row.names = FALSE, quote = TRUE)

original <- read.csv(original_module_path, check.names = FALSE)
original$gene <- toupper(as.character(original$gene))
shared <- merge(module_table, original[, c("gene", "module_color")], by = "gene")
names(shared)[names(shared) == "module_color"] <- "original42_module"

comb2 <- function(x) x * (x - 1) / 2
adjusted_rand_index <- function(tab) {
  tab <- as.matrix(tab)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  sum_ij <- sum(comb2(tab))
  sum_i <- sum(comb2(rowSums(tab)))
  sum_j <- sum(comb2(colSums(tab)))
  total <- comb2(n)
  expected <- sum_i * sum_j / total
  max_index <- (sum_i + sum_j) / 2
  if (max_index == expected) return(NA_real_)
  (sum_ij - expected) / (max_index - expected)
}
ari <- adjusted_rand_index(table(shared$baseline_module, shared$original42_module))

best_pairs <- do.call(rbind, lapply(sort(unique(shared$baseline_module)), function(mod) {
  genes_mod <- shared$gene[shared$baseline_module == mod]
  candidates <- lapply(sort(unique(shared$original42_module)), function(orig) {
    genes_orig <- shared$gene[shared$original42_module == orig]
    overlap <- length(intersect(genes_mod, genes_orig))
    union_n <- length(union(genes_mod, genes_orig))
    data.frame(
      baseline_module = mod,
      original42_module = orig,
      overlap_n = overlap,
      baseline_module_n = length(genes_mod),
      original42_module_n = length(genes_orig),
      jaccard = ifelse(union_n > 0, overlap / union_n, NA_real_)
    )
  })
  candidates <- do.call(rbind, candidates)
  candidates[order(-candidates$jaccard, -candidates$overlap_n, candidates$original42_module), ][1, ]
}))
write.csv(best_pairs, outputs[4], row.names = FALSE, quote = TRUE)

candidate10 <- c("CDK6", "PARP1", "KIT", "SYK", "HIF1A", "TOP2A", "CA2", "CD38", "TERT", "GSK3B")
candidate_mapping <- merge(data.frame(gene = candidate10), module_table, by = "gene", all.x = TRUE, sort = FALSE)
candidate_mapping <- merge(candidate_mapping, original[, c("gene", "module_color")], by = "gene", all.x = TRUE, sort = FALSE)
names(candidate_mapping)[names(candidate_mapping) == "module_color"] <- "original42_module"
candidate_mapping$module_match_exact <- candidate_mapping$baseline_module == candidate_mapping$original42_module
write.csv(candidate_mapping, outputs[5], row.names = FALSE, quote = TRUE)

elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
summary <- data.frame(
  metric = c(
    "n_profiles", "n_unique_participants", "n_HD", "n_SAA_baseline", "n_selected_genes",
    "n_baseline_modules_including_grey", "fixed_power", "reestimated_powerEstimate",
    "shared_genes_with_original42", "adjusted_rand_index_vs_original42",
    "candidate10_exact_module_match_n", "elapsed_seconds"
  ),
  value = c(
    nrow(datExpr), length(unique(baseline_meta$subject)), sum(baseline_meta$disease == "HD"),
    sum(baseline_meta$disease == "SAA"), ncol(datExpr), length(unique(module_colors)), 18,
    ifelse(length(sft$powerEstimate) == 0 || is.na(sft$powerEstimate), NA, sft$powerEstimate),
    nrow(shared), ari, sum(candidate_mapping$module_match_exact, na.rm = TRUE), elapsed
  )
)
write.table(summary, outputs[6], sep = "\t", row.names = FALSE, quote = FALSE)

report <- c(
  "# Baseline-only / one-profile-per-participant WGCNA sensitivity",
  "",
  paste0("Generated: ", Sys.time()),
  "",
  "## Design",
  "",
  "- 23 profiles, exactly one per participant: 19 SAA pretreatment profiles and 4 healthy donors.",
  "- The same collapsed counts, expression filter (detected in at least 8 profiles; total count at least 50; mitochondrial genes excluded), top-5,000 MAD selection, signed bicor network, fixed power 18, minimum module size 30, and merge cut height 0.25 were used.",
  "- This analysis removes repeated longitudinal profiles from network construction. It does not increase the number of independent healthy donors.",
  "",
  "## Results",
  "",
  paste0("- Baseline-only modules including grey: ", length(unique(module_colors)), "."),
  paste0("- Shared genes with the original 42-profile donor-timepoint network: ", nrow(shared), "."),
  paste0("- Adjusted Rand index of module assignments versus the original 42-profile network: ", format(ari, digits = 5), "."),
  paste0("- Frozen ten candidates with exactly matching color labels: ", sum(candidate_mapping$module_match_exact, na.rm = TRUE), "/10. Exact color-label agreement is a conservative diagnostic because WGCNA color names are nominal; the best-pair Jaccard table is the preferred comparison."),
  "",
  "## Interpretation boundary",
  "",
  "This run establishes computational feasibility and quantifies module reproducibility after retaining one profile per participant. With only four healthy donors, it is underpowered for strong between-group module-trait inference and should be reported as a sensitivity analysis, not used to replace the donor-aware 42-profile network or to claim independent validation.",
  "",
  paste0("Runtime: ", format(elapsed, digits = 6), " seconds.")
)
writeLines(report, outputs[7], useBytes = TRUE)
writeLines(c(
  paste0("started=", started),
  paste0("finished=", Sys.time()),
  paste0("elapsed_seconds=", elapsed),
  paste0("seed=20260902"),
  paste0("meta_input=", meta_path),
  paste0("counts_input=", counts_path),
  paste0("original_module_input=", original_module_path)
), outputs[8], useBytes = TRUE)

print(summary)
