# GSE165870 bulk differential-expression workflow
# Purpose:
#   Compare AA and healthy Lin-CD34+ HSPC bulk RNA-seq count profiles.
#   Intermediate tables and figures are retained to support reproducibility and quality control.
#
# 本脚本不做 WGCNA 主分析。
# 原因：
#   GSE165870 只有 3 个 healthy + 6 个 AA 样本，适合做 DEG / GSEA / 候选靶点表达验证，
#   但样本量偏小，不适合作为本文主 WGCNA 结论来源。

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
  library(readr)
})

# 1. 设置输入和输出路径 -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript gse165870_bulk_differential_expression.R <counts.txt> <candidate126.csv> <output_dir>")
}
count_file <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
candidate_file <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "00_GSE165870_bulk_DEG_run_log.txt")
sink(log_file, split = TRUE)
cat("GSE165870 bulk DEG analysis started at:", as.character(Sys.time()), "\n\n")

# 2. 读取 count 矩阵 ----------------------------------------------------------
# 原始矩阵结构：
#   GeneName | GeneSymbol | Ctrl2 | Ctrl7 | Ctrl8 | P18 | P19 | P20 | P21 | P22 | P23
# 其中 Ctrl 为 healthy，P 为 AA patient。
raw_counts <- read.delim(count_file, check.names = FALSE, stringsAsFactors = FALSE)
cat("Raw count matrix dimensions:", paste(dim(raw_counts), collapse = " x "), "\n")
cat("Raw columns:\n")
print(colnames(raw_counts))

# 3. 基因名整理与重复 GeneSymbol 合并 ----------------------------------------
# 同一个 GeneSymbol 可能对应多个 Ensembl ID。为了便于和候选靶点映射，这里按 GeneSymbol 合并 counts。
# 合并方式：同一 GeneSymbol 的 counts 求和。
sample_cols <- setdiff(colnames(raw_counts), c("GeneName", "GeneSymbol"))

counts_by_symbol <- raw_counts %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  group_by(GeneSymbol) %>%
  summarise(across(all_of(sample_cols), ~sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

count_mat <- as.data.frame(counts_by_symbol)
rownames(count_mat) <- count_mat$GeneSymbol
count_mat$GeneSymbol <- NULL
count_mat <- as.matrix(count_mat)
storage.mode(count_mat) <- "integer"

cat("Deduplicated count matrix dimensions:", paste(dim(count_mat), collapse = " x "), "\n")
write.csv(count_mat, file.path(out_dir, "01_GSE165870_counts_by_GeneSymbol.csv"), quote = FALSE)

# 4. 构建样本分组信息 ---------------------------------------------------------
# Ctrl2/Ctrl7/Ctrl8 为 healthy control；P18-P23 为 AA。
sample_id <- colnames(count_mat)
condition <- ifelse(grepl("^Ctrl", sample_id), "Healthy", "AA")
metadata <- data.frame(
  sample_id = sample_id,
  condition = factor(condition, levels = c("Healthy", "AA")),
  row.names = sample_id,
  stringsAsFactors = FALSE
)

cat("\nSample metadata:\n")
print(metadata)
write.csv(metadata, file.path(out_dir, "02_GSE165870_sample_metadata.csv"), quote = FALSE)

if (!identical(rownames(metadata), colnames(count_mat))) {
  stop("Sample order mismatch between metadata and count matrix.")
}

# 5. 低表达过滤和 DESeq2 差异分析 -------------------------------------------
# 过滤原则：保留至少在所有样本总 counts > 10 的基因，避免极低表达基因影响统计。
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = metadata,
  design = ~ condition
)

keep <- rowSums(counts(dds)) > 10
dds <- dds[keep, ]
cat("\nGenes retained after low-count filtering:", nrow(dds), "\n")

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "AA", "Healthy"))
res_df <- as.data.frame(res)
res_df$GeneSymbol <- rownames(res_df)
res_df <- res_df[, c("GeneSymbol", setdiff(colnames(res_df), "GeneSymbol"))]

# 标记差异方向：padj < 0.05 且 |log2FC| >= 1。
res_df$change <- "NOT"
res_df$change[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >= 1] <- "UP"
res_df$change[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange <= -1] <- "DOWN"

write.csv(res_df, file.path(out_dir, "03_DESeq2_all_genes_AA_vs_Healthy.csv"), row.names = FALSE, quote = FALSE)
sig_df <- res_df %>% filter(change %in% c("UP", "DOWN"))
write.csv(sig_df, file.path(out_dir, "04_DESeq2_significant_DEGs_padj0.05_log2FC1.csv"), row.names = FALSE, quote = FALSE)

cat("\nDEG summary:\n")
print(table(res_df$change))

# 6. 归一化表达矩阵和 VST 矩阵 ------------------------------------------------
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file.path(out_dir, "05_DESeq2_normalized_counts.csv"), quote = FALSE)

vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
write.csv(vst_mat, file.path(out_dir, "06_VST_expression_matrix.csv"), quote = FALSE)

# 7. PCA 图：检查 Healthy 和 AA 是否有整体表达分离 ---------------------------
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(size = 3.5, max.overlaps = 20) +
  scale_color_manual(values = c("Healthy" = "#4C78A8", "AA" = "#F58518")) +
  labs(
    title = "GSE165870 PCA: AA vs Healthy HSPCs",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(out_dir, "07_PCA_AA_vs_Healthy.png"), p_pca, width = 6.5, height = 5.2, dpi = 300)
ggsave(file.path(out_dir, "07_PCA_AA_vs_Healthy.pdf"), p_pca, width = 6.5, height = 5.2)

# 8. 样本相关性热图：检查样本间相似性 ----------------------------------------
sample_cor <- cor(vst_mat, method = "pearson")
write.csv(sample_cor, file.path(out_dir, "08_sample_correlation_matrix.csv"), quote = FALSE)

ann_col <- data.frame(condition = metadata$condition)
rownames(ann_col) <- rownames(metadata)

png(file.path(out_dir, "09_sample_correlation_heatmap.png"), width = 1800, height = 1500, res = 220)
pheatmap(
  sample_cor,
  annotation_col = ann_col,
  annotation_row = ann_col,
  color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
  main = "Sample correlation based on VST expression"
)
dev.off()

pdf(file.path(out_dir, "09_sample_correlation_heatmap.pdf"), width = 7, height = 6)
pheatmap(
  sample_cor,
  annotation_col = ann_col,
  annotation_row = ann_col,
  color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
  main = "Sample correlation based on VST expression"
)
dev.off()

# 9. 火山图：展示 AA vs Healthy 差异基因 -------------------------------------
volcano_df <- res_df %>%
  mutate(
    neg_log10_padj = -log10(padj),
    neg_log10_padj = ifelse(is.infinite(neg_log10_padj), NA, neg_log10_padj)
  )

top_labels <- volcano_df %>%
  filter(change %in% c("UP", "DOWN")) %>%
  arrange(padj) %>%
  head(12)

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj, color = change)) +
  geom_point(alpha = 0.65, size = 1.6, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.35) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.35) +
  geom_text_repel(data = top_labels, aes(label = GeneSymbol), size = 3, max.overlaps = 30) +
  scale_color_manual(values = c("DOWN" = "#4C78A8", "NOT" = "grey78", "UP" = "#F58518")) +
  labs(
    title = "Differential expression in GSE165870",
    x = "log2 fold change (AA / Healthy)",
    y = "-log10 adjusted P value"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(out_dir, "10_volcano_AA_vs_Healthy.png"), p_volcano, width = 6.5, height = 5.6, dpi = 300)
ggsave(file.path(out_dir, "10_volcano_AA_vs_Healthy.pdf"), p_volcano, width = 6.5, height = 5.6)

# 10. Top DEG 热图 ------------------------------------------------------------
# 选择 padj 最小的前 50 个 DEG。如果显著 DEG 少于 50，则取全部显著 DEG。
top_heatmap_genes <- sig_df %>%
  arrange(padj) %>%
  head(50) %>%
  pull(GeneSymbol)

if (length(top_heatmap_genes) >= 2) {
  top_mat <- vst_mat[top_heatmap_genes, , drop = FALSE]
  top_mat_z <- t(scale(t(top_mat)))
  top_mat_z[is.na(top_mat_z)] <- 0
  write.csv(top_mat_z, file.path(out_dir, "11_top_DEG_heatmap_zscore_matrix.csv"), quote = FALSE)
  
  png(file.path(out_dir, "12_top_DEG_heatmap.png"), width = 1600, height = 2000, res = 220)
  pheatmap(
    top_mat_z,
    annotation_col = ann_col,
    show_rownames = TRUE,
    fontsize_row = 6,
    color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
    main = "Top DEGs in GSE165870"
  )
  dev.off()
  
  pdf(file.path(out_dir, "12_top_DEG_heatmap.pdf"), width = 6.5, height = 8)
  pheatmap(
    top_mat_z,
    annotation_col = ann_col,
    show_rownames = TRUE,
    fontsize_row = 6,
    color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
    main = "Top DEGs in GSE165870"
  )
  dev.off()
} else {
  cat("\nTop DEG heatmap skipped because fewer than 2 significant DEGs were detected.\n")
}

# 11. Map the prespecified 126-candidate set to GSE165870 -------------------
# Differential expression provides complementary HSPC bulk-expression evidence
# and does not define or rerank the candidate set.
candidate_targets <- read.csv(candidate_file, stringsAsFactors = FALSE)
candidate_genes <- unique(candidate_targets$GeneSymbol)

candidate_map <- res_df %>%
  filter(GeneSymbol %in% candidate_genes) %>%
  arrange(padj)

candidate_map$expressed_in_GSE165870 <- candidate_map$baseMean > 0
candidate_map$is_DEG_padj0.05_log2FC1 <- candidate_map$change %in% c("UP", "DOWN")

missing_candidates <- setdiff(candidate_genes, candidate_map$GeneSymbol)
missing_df <- data.frame(
  GeneSymbol = missing_candidates,
  baseMean = NA,
  log2FoldChange = NA,
  lfcSE = NA,
  stat = NA,
  pvalue = NA,
  padj = NA,
  change = "NOT_IN_MATRIX",
  expressed_in_GSE165870 = FALSE,
  is_DEG_padj0.05_log2FC1 = FALSE
)

candidate_map_full <- bind_rows(candidate_map, missing_df) %>%
  arrange(desc(expressed_in_GSE165870), padj)

write.csv(candidate_map_full, file.path(out_dir, "13_candidate126_mapped_to_GSE165870_DESeq2.csv"), row.names = FALSE, quote = FALSE)

cat("\nCandidate target mapping summary:\n")
cat("Input candidate targets:", length(candidate_genes), "\n")
cat("Found in DESeq2 matrix:", nrow(candidate_map), "\n")
cat("Missing from DESeq2 matrix:", length(missing_candidates), "\n")
cat("Candidate targets meeting DEG cutoff:", sum(candidate_map_full$is_DEG_padj0.05_log2FC1, na.rm = TRUE), "\n")

# 12. 候选靶点表达热图 --------------------------------------------------------
# 展示 126 候选靶点中，在本数据集表达且变化最明显的前 40 个。
candidate_heatmap_genes <- candidate_map %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(40) %>%
  pull(GeneSymbol)

if (length(candidate_heatmap_genes) >= 2) {
  cand_mat <- vst_mat[candidate_heatmap_genes, , drop = FALSE]
  cand_mat_z <- t(scale(t(cand_mat)))
  cand_mat_z[is.na(cand_mat_z)] <- 0
  write.csv(cand_mat_z, file.path(out_dir, "14_candidate_target_heatmap_zscore_matrix.csv"), quote = FALSE)
  
  png(file.path(out_dir, "15_candidate_target_heatmap.png"), width = 1600, height = 1900, res = 220)
  pheatmap(
    cand_mat_z,
    annotation_col = ann_col,
    show_rownames = TRUE,
    fontsize_row = 7,
    color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
    main = "Candidate targets mapped to GSE165870"
  )
  dev.off()
  
  pdf(file.path(out_dir, "15_candidate_target_heatmap.pdf"), width = 6.5, height = 7.5)
  pheatmap(
    cand_mat_z,
    annotation_col = ann_col,
    show_rownames = TRUE,
    fontsize_row = 7,
    color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
    main = "Candidate targets mapped to GSE165870"
  )
  dev.off()
}

# 13. Save analysis objects for reproducibility -----------------------------
save(dds, res_df, sig_df, norm_counts, vst_mat, metadata, candidate_map_full,
     file = file.path(out_dir, "16_GSE165870_DESeq2_analysis_objects.RData"))

cat("\nAnalysis finished at:", as.character(Sys.time()), "\n")
sink()
