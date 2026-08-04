suppressPackageStartupMessages({
  library(WGCNA)
  library(limma)
})

options(stringsAsFactors = FALSE)
set.seed(20260715)
disableWGCNAThreads()
allowWGCNAThreads(nThreads = 4)

bundle_root <- normalizePath(if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[1] else ".", mustWork = TRUE)
out_dir <- file.path(bundle_root, "data", "wgcna_strict", "outputs")
input_dir <- file.path(bundle_root, "data", "wgcna_strict", "inputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- "strict_donor_aware_wgcna"

input_files <- c(
  metadata = file.path(input_dir, "original_sample_metadata.tsv"),
  counts = file.path(input_dir, "original_pseudobulk_counts.csv"),
  original_module_gene_kME_GS = file.path(input_dir, "original_module_gene_table.csv"),
  original_candidate_mapping = file.path(input_dir, "frozen_candidate_pool_mapping.csv"),
  original_top30_scores = file.path(input_dir, "frozen_top30_reference_scores.csv"),
  candidate_pool = file.path(input_dir, "candidate_pool_126.txt")
)

planned_outputs <- file.path(out_dir, c(
  paste0(prefix, "_run_log.txt"),
  paste0(prefix, "_input_sha256.tsv"),
  paste0(prefix, "_collapsed_sample_metadata.tsv"),
  paste0(prefix, "_collapsed_counts_by_subject_timepoint.csv"),
  paste0(prefix, "_collapsed_log1pCP10K_by_subject_timepoint.csv"),
  paste0(prefix, "_gene_filtering_record.csv"),
  paste0(prefix, "_soft_threshold_fit.csv"),
  paste0(prefix, "_module_gene_table.csv"),
  paste0(prefix, "_module_signature_overlap_annotation.csv"),
  paste0(prefix, "_strict_group_contrast_module_tests.csv"),
  paste0(prefix, "_candidate126_strict_wgcna_evidence_long.csv"),
  paste0(prefix, "_candidate126_strict_wgcna_evidence_summary.csv"),
  paste0(prefix, "_original_top30_strict_status.csv"),
  paste0(prefix, "_original48_vs_42_module_overlap_matrix.csv"),
  paste0(prefix, "_original48_vs_42_module_overlap_best_pairs.csv"),
  paste0(prefix, "_module_overlap_ARI.txt"),
  paste0(prefix, "_session_info.txt"),
  paste0(prefix, "_QC_report.md"),
  paste0(prefix, "_final_all_outputs_sha256_20260715.tsv")
))
existing_outputs <- planned_outputs[file.exists(planned_outputs)]
if (length(existing_outputs) > 0) {
  stop("Refusing to overwrite existing 92 outputs: ", paste(existing_outputs, collapse = "; "))
}

log_file <- file.path(out_dir, paste0(prefix, "_run_log.txt"))
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("CD34 strict donor-aware WGCNA sensitivity started:", as.character(Sys.time()), "\n")
cat("Prefix:", prefix, "\n")
cat("Seed: 20260715\n\n")

write_csv <- function(x, filename, row.names = FALSE) {
  write.csv(x, file.path(out_dir, filename), row.names = row.names, quote = TRUE, fileEncoding = "UTF-8")
}
write_tsv <- function(x, filename) {
  write.table(x, file.path(out_dir, filename), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
}
read_gene_matrix <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  gene <- toupper(x[[1]])
  mat <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- gene
  if (any(duplicated(rownames(mat)))) {
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  }
  mat
}
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], use = "pairwise.complete.obs"))
}
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

signature_sets <- list(
  HSPC_identity = c("CD34", "KIT", "PROM1", "GATA2", "MECOM", "RUNX1", "HLF", "MEIS1"),
  HSPC_injury_stress = c("TP53", "CDKN1A", "BAX", "BCL2L11", "CASP3", "FAS", "DDIT3", "ATF4", "GADD45A", "BBC3"),
  immune_IFN_TNF = c("STAT1", "IRF1", "ISG15", "IFIT1", "IFIT3", "CXCL10", "TNF", "NFKB1", "RELA", "JUN", "FOS"),
  hematopoietic_support = c("MPL", "JAK2", "STAT5A", "STAT5B", "PIM1", "BCL2L1", "MYC", "IGF1R", "MET", "CXCR4", "KIT"),
  liver_bone_marrow_response = c("MPL", "JAK2", "STAT5A", "STAT5B", "IGF1R", "MET", "HAMP", "FTH1", "FTL", "TFRC", "SLC40A1"),
  cell_cycle_recovery = c("MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "CDK1", "CCNB1", "TYMS")
)

input_sha <- data.frame(
  label = names(input_files),
  path = unname(input_files),
  exists = file.exists(input_files),
  sha256 = ifelse(file.exists(input_files), unname(tools::sha256sum(input_files)), NA_character_),
  stringsAsFactors = FALSE
)
write_tsv(input_sha, paste0(prefix, "_input_sha256.tsv"))
print(input_sha)

meta <- read.delim(input_files[["metadata"]], stringsAsFactors = FALSE, check.names = FALSE)
counts <- read_gene_matrix(input_files[["counts"]])
candidate_genes <- unique(toupper(trimws(readLines(input_files[["candidate_pool"]], warn = FALSE, encoding = "UTF-8"))))
candidate_genes <- candidate_genes[nzchar(candidate_genes)]
original_candidate <- read.csv(input_files[["original_candidate_mapping"]], stringsAsFactors = FALSE, check.names = FALSE)
original_top30 <- read.csv(input_files[["original_top30_scores"]], stringsAsFactors = FALSE, check.names = FALSE)
original_gene_table <- read.csv(input_files[["original_module_gene_kME_GS"]], stringsAsFactors = FALSE, check.names = FALSE)

meta$count_suffix <- as.character(meta$count_suffix)
meta$subject <- as.character(meta$subject)
meta$timepoint_collapsed <- ifelse(meta$disease == "HD" | !nzchar(meta$timepoint), "HD", meta$timepoint)
meta$subject_timepoint_id <- paste(meta$subject, meta$timepoint_collapsed, sep = "__")
stopifnot(all(meta$count_suffix %in% colnames(counts)))
counts <- counts[, meta$count_suffix, drop = FALSE]

key_levels <- unique(meta$subject_timepoint_id)
collapsed_counts <- sapply(key_levels, function(k) {
  suffixes <- meta$count_suffix[meta$subject_timepoint_id == k]
  rowSums(counts[, suffixes, drop = FALSE])
})
collapsed_counts <- as.matrix(collapsed_counts)
storage.mode(collapsed_counts) <- "numeric"
rownames(collapsed_counts) <- rownames(counts)

collapsed_meta <- do.call(rbind, lapply(key_levels, function(k) {
  sub <- meta[meta$subject_timepoint_id == k, , drop = FALSE]
  matrix_umi <- colSums(counts[, sub$count_suffix, drop = FALSE])
  data.frame(
    subject_timepoint_id = k,
    subject = sub$subject[1],
    disease = sub$disease[1],
    timepoint = sub$timepoint_collapsed[1],
    group = ifelse(sub$disease[1] == "HD", "HD", paste0("SAA_", sub$timepoint_collapsed[1])),
    n_input_profiles = nrow(sub),
    input_count_suffixes = paste(sub$count_suffix, collapse = ";"),
    input_geo_accessions = paste(sub$geo_accession, collapse = ";"),
    had_technical_repeat = nrow(sub) > 1 || any(tolower(as.character(sub$is_repeat)) %in% c("true", "1", "yes")),
    n_repeat_flagged_profiles = sum(tolower(as.character(sub$is_repeat)) %in% c("true", "1", "yes")),
    n_cells_sum = sum(as.numeric(sub$n_cells), na.rm = TRUE),
    total_umi_metadata_sum = sum(as.numeric(sub$total_umi), na.rm = TRUE),
    total_umi_from_matrix_sum = sum(matrix_umi, na.rm = TRUE),
    total_umi_delta_sum = sum(matrix_umi, na.rm = TRUE) - sum(as.numeric(sub$total_umi), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write_tsv(collapsed_meta, paste0(prefix, "_collapsed_sample_metadata.tsv"))
write_csv(data.frame(gene = rownames(collapsed_counts), collapsed_counts, check.names = FALSE), paste0(prefix, "_collapsed_counts_by_subject_timepoint.csv"))

cat("Original rows:", nrow(meta), "\n")
cat("Collapsed subject x timepoint rows:", nrow(collapsed_meta), "\n")
cat("Subjects:", length(unique(collapsed_meta$subject)), "\n")
print(table(collapsed_meta$group))

lib_size <- colSums(collapsed_counts)
log1p_cp10k <- log1p(sweep(collapsed_counts, 2, lib_size / 10000, "/"))
log1p_cp10k[!is.finite(log1p_cp10k)] <- 0
write_csv(data.frame(gene = rownames(log1p_cp10k), log1p_cp10k, check.names = FALSE), paste0(prefix, "_collapsed_log1pCP10K_by_subject_timepoint.csv"))

detected_samples <- rowSums(collapsed_counts > 0)
total_counts <- rowSums(collapsed_counts)
is_mito <- grepl("^MT-", rownames(collapsed_counts))
min_detected_samples <- 8
keep_expr <- detected_samples >= min_detected_samples & total_counts >= 50 & !is_mito
mad_value <- apply(log1p_cp10k[keep_expr, , drop = FALSE], 1, mad, na.rm = TRUE)
mad_value[!is.finite(mad_value)] <- 0
wgcna_genes <- names(sort(mad_value, decreasing = TRUE))[seq_len(min(5000, length(mad_value)))]

filter_table <- data.frame(
  gene = rownames(collapsed_counts),
  detected_subject_timepoints = detected_samples,
  total_counts = total_counts,
  is_mito = is_mito,
  pass_expression_filter = keep_expr,
  mad_log1p_cp10k = ifelse(rownames(collapsed_counts) %in% names(mad_value), mad_value[rownames(collapsed_counts)], NA_real_),
  selected_for_WGCNA = rownames(collapsed_counts) %in% wgcna_genes,
  fixed_min_detected_subject_timepoints = min_detected_samples,
  stringsAsFactors = FALSE
)
write_csv(filter_table, paste0(prefix, "_gene_filtering_record.csv"))

datExpr0 <- as.data.frame(t(log1p_cp10k[wgcna_genes, collapsed_meta$subject_timepoint_id, drop = FALSE]), check.names = FALSE)
rownames(datExpr0) <- collapsed_meta$subject_timepoint_id
subject_block <- collapsed_meta$subject[match(rownames(datExpr0), collapsed_meta$subject_timepoint_id)]
group <- factor(collapsed_meta$group[match(rownames(datExpr0), collapsed_meta$subject_timepoint_id)], levels = c("HD", "SAA_baseline", "SAA_3M", "SAA_6M"))

gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  subject_block <- subject_block[gsg$goodSamples]
  group <- group[gsg$goodSamples]
}

powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(
  datExpr0,
  powerVector = powers,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = list(maxPOutliers = 0.1, use = "pairwise.complete.obs"),
  verbose = 5
)
fit_indices <- sft$fitIndices
fit_indices$strict_fixed_power_used <- 18
fit_indices$reestimated_powerEstimate <- sft$powerEstimate
write_csv(fit_indices, paste0(prefix, "_soft_threshold_fit.csv"))
cat("Fixed WGCNA power used: 18\n")
cat("pickSoftThreshold powerEstimate:", sft$powerEstimate, "\n")

net <- blockwiseModules(
  datExpr0,
  power = 18,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  maxPOutliers = 0.1,
  maxBlockSize = ncol(datExpr0),
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)
moduleColors <- net$colors
names(moduleColors) <- colnames(datExpr0)
MEs <- orderMEs(net$MEs)
rownames(MEs) <- rownames(datExpr0)
kME <- signedKME(datExpr0, MEs, outputColumnName = "kME")

module_gene_table <- data.frame(
  gene = names(moduleColors),
  module_color = unname(moduleColors),
  stringsAsFactors = FALSE
)
for (mod in sort(unique(moduleColors))) {
  col <- paste0("kME", mod)
  module_gene_table[[col]] <- if (col %in% colnames(kME)) kME[module_gene_table$gene, col] else NA_real_
}
write_csv(module_gene_table, paste0(prefix, "_module_gene_table.csv"))

module_signature_overlap <- do.call(rbind, lapply(sort(unique(moduleColors)), function(mod) {
  genes <- names(moduleColors)[moduleColors == mod]
  overlaps <- sapply(signature_sets, function(sig) length(intersect(genes, toupper(sig))))
  data.frame(
    module_color = mod,
    module_eigengene = paste0("ME", mod),
    n_genes = length(genes),
    t(as.data.frame(overlaps)),
    candidate126_genes_in_module = paste(intersect(genes, candidate_genes), collapse = ";"),
    n_candidate126_in_module = length(intersect(genes, candidate_genes)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))
write_csv(module_signature_overlap, paste0(prefix, "_module_signature_overlap_annotation.csv"))

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
expr_modules <- t(as.matrix(MEs))
dup_cor <- duplicateCorrelation(expr_modules, design, block = subject_block)
fit <- lmFit(expr_modules, design, block = subject_block, correlation = dup_cor$consensus.correlation)
contrast_matrix <- makeContrasts(
  baseline_minus_HD = SAA_baseline - HD,
  M3_minus_baseline = SAA_3M - SAA_baseline,
  M6_minus_baseline = SAA_6M - SAA_baseline,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

contrast_tests <- do.call(rbind, lapply(colnames(contrast_matrix), function(contrast_name) {
  tt <- topTable(fit2, coef = contrast_name, number = Inf, sort.by = "none")
  tt$module <- rownames(tt)
  tt$contrast <- contrast_name
  tt$module_color <- sub("^ME", "", tt$module)
  tt$duplicate_correlation <- dup_cor$consensus.correlation
  tt$n_subject_timepoints <- nrow(MEs)
  tt$n_subjects <- length(unique(subject_block))
  tt$raw_module_group_contrast_effect <- tt$logFC
  tt[, c(
    "module", "module_color", "contrast", "logFC", "t", "P.Value", "adj.P.Val", "B",
    "raw_module_group_contrast_effect", "duplicate_correlation", "n_subject_timepoints", "n_subjects"
  )]
}))
contrast_tests$BH_by_contrast <- ave(contrast_tests$P.Value, contrast_tests$contrast, FUN = function(x) p.adjust(x, method = "BH"))
contrast_tests$BH_global_all_contrasts <- p.adjust(contrast_tests$P.Value, method = "BH")
contrast_tests <- contrast_tests[order(contrast_tests$contrast, contrast_tests$P.Value), ]
write_csv(contrast_tests, paste0(prefix, "_strict_group_contrast_module_tests.csv"))

module_contrast_lookup <- contrast_tests
candidate_rows <- do.call(rbind, lapply(candidate_genes, function(gene) {
  in_expr <- gene %in% rownames(log1p_cp10k)
  in_top5000 <- gene %in% names(moduleColors)
  own_module <- if (in_top5000) unname(moduleColors[gene]) else NA_character_
  own_kme <- NA_real_
  if (in_top5000) {
    kme_col <- paste0("kME", own_module)
    if (kme_col %in% colnames(kME)) own_kme <- unname(kME[gene, kme_col])
  }
  projection_module <- NA_character_
  projection_cor <- NA_real_
  if (in_expr) {
    expr <- as.numeric(log1p_cp10k[gene, rownames(MEs)])
    cors <- apply(MEs, 2, safe_cor, y = expr)
    if (any(is.finite(cors))) {
      best <- names(cors)[which.max(abs(cors))]
      projection_module <- sub("^ME", "", best)
      projection_cor <- unname(cors[best])
    }
  }
  module_for_context <- if (in_top5000) own_module else projection_module
  ct <- module_contrast_lookup[module_contrast_lookup$module_color == module_for_context, , drop = FALSE]
  if (nrow(ct) == 0) {
    ct <- data.frame(
      contrast = colnames(contrast_matrix),
      logFC = NA_real_,
      t = NA_real_,
      P.Value = NA_real_,
      BH_by_contrast = NA_real_,
      BH_global_all_contrasts = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, lapply(seq_len(nrow(ct)), function(i) {
    data.frame(
      gene = gene,
      expressed_in_collapsed_matrix = in_expr,
      eligible_for_strict_wgcna_evidence = in_top5000,
      evidence_reason = ifelse(in_top5000, "Top5000 gene: own module kME plus prespecified group contrast statistics", "Not Top5000: projection context only, excluded from strict WGCNA evidence rank"),
      own_module_if_top5000 = own_module,
      abs_own_kME_if_top5000 = abs(own_kme),
      projection_context_module = projection_module,
      projection_context_cor = projection_cor,
      contrast = ct$contrast[i],
      module_contrast_logFC = ct$logFC[i],
      module_contrast_t = ct$t[i],
      module_contrast_p = ct$P.Value[i],
      module_contrast_BH_by_contrast = ct$BH_by_contrast[i],
      module_contrast_BH_global = ct$BH_global_all_contrasts[i],
      stringsAsFactors = FALSE
    )
  }))
}))
candidate_rows$strict_support_level <- ifelse(
  candidate_rows$eligible_for_strict_wgcna_evidence &
    is.finite(candidate_rows$abs_own_kME_if_top5000) &
    candidate_rows$abs_own_kME_if_top5000 >= 0.35 &
    is.finite(candidate_rows$module_contrast_BH_by_contrast) &
    candidate_rows$module_contrast_BH_by_contrast < 0.05,
  "FDR_support",
  ifelse(
    candidate_rows$eligible_for_strict_wgcna_evidence &
      is.finite(candidate_rows$abs_own_kME_if_top5000) &
      candidate_rows$abs_own_kME_if_top5000 >= 0.35 &
      is.finite(candidate_rows$module_contrast_p) &
      candidate_rows$module_contrast_p < 0.05,
    "nominal_only",
    ifelse(candidate_rows$eligible_for_strict_wgcna_evidence, "no_support", "projection_only_not_ranked")
  )
)
write_csv(candidate_rows, paste0(prefix, "_candidate126_strict_wgcna_evidence_long.csv"))

candidate_summary <- do.call(rbind, lapply(split(candidate_rows, candidate_rows$gene), function(df) {
  support_order <- c(FDR_support = 1, nominal_only = 2, no_support = 3, projection_only_not_ranked = 4)
  df$ord <- support_order[df$strict_support_level]
  best <- df[order(df$ord, df$module_contrast_BH_by_contrast, df$module_contrast_p, -abs(df$module_contrast_t)), ][1, ]
  data.frame(
    gene = best$gene,
    expressed_in_collapsed_matrix = best$expressed_in_collapsed_matrix,
    eligible_for_strict_wgcna_evidence = best$eligible_for_strict_wgcna_evidence,
    own_module_if_top5000 = best$own_module_if_top5000,
    abs_own_kME_if_top5000 = best$abs_own_kME_if_top5000,
    projection_context_module = best$projection_context_module,
    projection_context_cor = best$projection_context_cor,
    best_strict_support_level = best$strict_support_level,
    best_prespecified_contrast = best$contrast,
    best_module_contrast_logFC = best$module_contrast_logFC,
    best_module_contrast_t = best$module_contrast_t,
    best_module_contrast_p = best$module_contrast_p,
    best_module_contrast_BH_by_contrast = best$module_contrast_BH_by_contrast,
    best_module_contrast_BH_global = best$module_contrast_BH_global,
    n_FDR_supported_contrasts = sum(df$strict_support_level == "FDR_support"),
    n_nominal_supported_contrasts = sum(df$strict_support_level %in% c("FDR_support", "nominal_only")),
    stringsAsFactors = FALSE
  )
}))
candidate_summary$strict_evidence_order <- NA_integer_
rankable <- candidate_summary$eligible_for_strict_wgcna_evidence
support_rank <- c(FDR_support = 1, nominal_only = 2, no_support = 3, projection_only_not_ranked = 4)
ord <- order(
  support_rank[candidate_summary$best_strict_support_level[rankable]],
  candidate_summary$best_module_contrast_BH_by_contrast[rankable],
  candidate_summary$best_module_contrast_p[rankable],
  -candidate_summary$abs_own_kME_if_top5000[rankable],
  candidate_summary$gene[rankable]
)
candidate_summary$strict_evidence_order[which(rankable)[ord]] <- seq_len(sum(rankable))
candidate_summary <- candidate_summary[order(is.na(candidate_summary$strict_evidence_order), candidate_summary$strict_evidence_order, candidate_summary$gene), ]
write_csv(candidate_summary, paste0(prefix, "_candidate126_strict_wgcna_evidence_summary.csv"))

top30_genes <- head(original_top30$GeneSymbol[order(suppressWarnings(as.numeric(original_top30$rank_network_integrated)))], 30)
top30 <- original_top30[match(top30_genes, original_top30$GeneSymbol), ]
top30_status <- merge(top30, candidate_summary, by.x = "GeneSymbol", by.y = "gene", all.x = TRUE, sort = FALSE)
top30_status$strict_top30_status <- top30_status$best_strict_support_level
top30_status$strict_status_interpretation <- ifelse(
  top30_status$strict_top30_status == "FDR_support",
  "FDR-supported module contrast; can be cited as strict donor-aware WGCNA support",
  ifelse(
    top30_status$strict_top30_status == "nominal_only",
    "Nominal module contrast only; use as exploratory WGCNA context, not robust support",
    ifelse(
      top30_status$strict_top30_status == "no_support",
      "Top5000 module membership but no prespecified contrast support",
      "Projection context only or absent from strict WGCNA evidence rank"
    )
  )
)
write_csv(top30_status, paste0(prefix, "_original_top30_strict_status.csv"))

original_modules <- original_gene_table[, c("gene", "module_color")]
original_modules$gene <- toupper(original_modules$gene)
names(original_modules)[2] <- "original48_module"
new_modules <- data.frame(gene = names(moduleColors), strict42_module = unname(moduleColors), stringsAsFactors = FALSE)
shared <- merge(original_modules, new_modules, by = "gene")
overlap_tab <- table(shared$original48_module, shared$strict42_module)
write_csv(as.data.frame.matrix(overlap_tab), paste0(prefix, "_original48_vs_42_module_overlap_matrix.csv"), row.names = TRUE)
ari <- adjusted_rand_index(overlap_tab)
best_pairs <- do.call(rbind, lapply(rownames(overlap_tab), function(orig) {
  counts <- overlap_tab[orig, ]
  best_new <- names(counts)[which.max(counts)]
  union_n <- sum(overlap_tab[orig, ]) + sum(overlap_tab[, best_new]) - max(counts)
  data.frame(
    original48_module = orig,
    best_strict42_module = best_new,
    overlap_genes = as.integer(max(counts)),
    original48_module_genes_in_shared = as.integer(sum(overlap_tab[orig, ])),
    strict42_module_genes_in_shared = as.integer(sum(overlap_tab[, best_new])),
    jaccard_on_shared_genes = ifelse(union_n > 0, as.numeric(max(counts)) / union_n, NA_real_),
    stringsAsFactors = FALSE
  )
}))
best_pairs <- best_pairs[order(-best_pairs$overlap_genes, -best_pairs$jaccard_on_shared_genes), ]
write_csv(best_pairs, paste0(prefix, "_original48_vs_42_module_overlap_best_pairs.csv"))
writeLines(c(
  paste("shared_genes", nrow(shared), sep = "\t"),
  paste("adjusted_rand_index", ari, sep = "\t"),
  "Note\tARI/overlap compares color partitions on shared Top5000 genes; color names are labels and should not be interpreted biologically."
), con = file.path(out_dir, paste0(prefix, "_module_overlap_ARI.txt")), useBytes = TRUE)

session_path <- file.path(out_dir, paste0(prefix, "_session_info.txt"))
writeLines(c(
  paste("Run time:", as.character(Sys.time())),
  paste("R:", R.version.string),
  paste("WGCNA:", as.character(packageVersion("WGCNA"))),
  paste("limma:", as.character(packageVersion("limma"))),
  paste("dynamicTreeCut:", as.character(packageVersion("dynamicTreeCut"))),
  paste("fastcluster:", as.character(packageVersion("fastcluster"))),
  paste("Seed:", 20260715),
  "",
  capture.output(sessionInfo())
), con = session_path, useBytes = TRUE)

status_counts <- table(top30_status$strict_top30_status, useNA = "ifany")
qc_lines <- c(
  "# QC1 strict donor-aware CD34 WGCNA sensitivity 20260715",
  "",
  "This participant-aware sensitivity analysis reports module-level evidence and does not alter the prespecified candidate order.",
  "",
  "## Strict inferential design",
  "",
  "- Only prespecified external clinical group contrasts were used for module inference: SAA_baseline - HD, SAA_3M - SAA_baseline, and SAA_6M - SAA_baseline.",
  "- Module tests used limma duplicateCorrelation/lmFit with subject block.",
  "- Signature scores and marker signatures were used only as overlap annotations, not as inferential traits, evidence scores, or post hoc minimum-P trait choices.",
  "- No artificial axis weights or minmax-scaled raw top-trait correlations were used for candidate evidence.",
  "",
  "## Independence and WGCNA",
  "",
  paste0("- Original profiles: ", nrow(meta), "."),
  paste0("- Collapsed subject x timepoint profiles: ", nrow(collapsed_meta), "."),
  paste0("- Subjects: ", length(unique(collapsed_meta$subject)), "."),
  paste0("- Group counts after collapse: ", paste(names(table(collapsed_meta$group)), as.integer(table(collapsed_meta$group)), sep = "=", collapse = "; "), "."),
  "- WGCNA used signed bicor network with fixed power=18, minModuleSize=30 and mergeCutHeight=0.25.",
  paste0("- pickSoftThreshold powerEstimate recorded as: ", sft$powerEstimate, "."),
  "",
  "## Top30 strict status",
  "",
  paste(capture.output(print(status_counts)), collapse = "\n"),
  "",
  "Interpretation:",
  "",
  "- FDR_support: Top5000 candidate whose own WGCNA module has BH<0.05 in at least one prespecified clinical contrast and |kME|>=0.35.",
  "- nominal_only: Top5000 candidate whose module has nominal P<0.05 but no contrast-level BH<0.05.",
  "- no_support: Top5000 module membership without prespecified contrast support.",
  "- projection_only_not_ranked: non-Top5000 candidate; projection is context only and excluded from strict WGCNA evidence rank.",
  "",
  "Only results from the participant-aware workflow are used for the reported module-level sensitivity analysis.",
  "",
  "## Module comparison",
  "",
  paste0("- Shared genes used for original48 vs strict42 module partition comparison: ", nrow(shared), "."),
  paste0("- Adjusted Rand index: ", signif(ari, 4), "."),
  "- Color changes are not biological interpretations; use the overlap matrix and ARI only for partition comparison.",
  "",
  "## Key outputs",
  "",
  paste0("- Module contrast tests: ", paste0(prefix, "_strict_group_contrast_module_tests.csv")),
  paste0("- Candidate strict evidence summary: ", paste0(prefix, "_candidate126_strict_wgcna_evidence_summary.csv")),
  paste0("- Original Top30 strict status: ", paste0(prefix, "_original_top30_strict_status.csv")),
  paste0("- Original48 vs strict42 overlap matrix: ", paste0(prefix, "_original48_vs_42_module_overlap_matrix.csv")),
  paste0("- Final SHA256 manifest: ", paste0(prefix, "_final_all_outputs_sha256_20260715.tsv")),
  "",
  "No original files were deleted, moved or overwritten."
)
writeLines(qc_lines, con = file.path(out_dir, paste0(prefix, "_QC_report.md")), useBytes = TRUE)

output_paths <- file.path(out_dir, list.files(out_dir, pattern = paste0("^", prefix), full.names = FALSE))
output_sha <- data.frame(
  file = basename(output_paths),
  length = file.info(output_paths)$size,
  sha256 = unname(tools::sha256sum(output_paths)),
  path = output_paths,
  stringsAsFactors = FALSE
)
write_tsv(output_sha, paste0(prefix, "_final_all_outputs_sha256_20260715.tsv"))

cat("\nCompleted strict donor-aware sensitivity:", as.character(Sys.time()), "\n")
cat("Generated prefix:", prefix, "\n")
