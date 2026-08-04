# GSE247531 BM official-style Seurat RPCA integration and split UMAP
#
# This script intentionally follows the official archived GSE247531 BM workflow:
#   NormalizeData -> FindVariableFeatures -> SelectIntegrationFeatures ->
#   ScaleData -> RunPCA -> FindIntegrationAnchors(reference = HD samples, RPCA) ->
#   IntegrateData -> RunPCA -> RunUMAP -> FindNeighbors -> FindClusters.
#
# Local adaptation:
#   1) GSE247531_BMcounts.Rdata.gz is rebuilt into per-sample Seurat objects.
#   2) Broad BM cell types are assigned by a transparent marker-panel score so
#      the UMAP can be split by time point/group and colored by cell compartment.
#   3) This is a visualization/context analysis for the BM microenvironment
#      story; it is not used alone to declare final core targets.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(flag) flag %in% args
get_arg <- function(prefix, default = "") {
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[1], fixed = TRUE)
}

dry_run <- has_flag("--dry-run")
project_root <- get_arg("--project-root=", ".")
out_dir <- get_arg("--out-dir=", file.path(project_root, "06_results/GSE247531_BM_official_style_Seurat_RPCA_UMAP"))
n_threads <- as.integer(get_arg("--threads=", "16"))
anchor_workers <- as.integer(get_arg("--anchor-workers=", as.character(min(n_threads, 4))))
integrate_workers <- as.integer(get_arg("--integrate-workers=", "1"))
dims_use <- seq_len(as.integer(get_arg("--dims=", "50")))
resolution_values <- as.numeric(strsplit(get_arg("--resolutions=", "0.8,0.5,1,1.5"), ",", fixed = TRUE)[[1]])
save_object <- tolower(get_arg("--save-object=", "true")) %in% c("true", "t", "1", "yes")

raw_dir <- get_arg("--raw-dir=", file.path(project_root, "external_data/GSE247531"))
bm_file <- file.path(raw_dir, "GSE247531_BMcounts.Rdata.gz")
expected_bm_file_size <- 1821662122
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "97_BM_official_style_Seurat_RPCA_UMAP_run_log.txt")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("BM official-style Seurat RPCA UMAP started at:", as.character(Sys.time()), "\n")
cat("Project root:", project_root, "\n")
cat("BM file:", bm_file, "\n")
cat("Output directory:", out_dir, "\n")
cat("Dry run:", dry_run, "\n")
cat("Threads:", n_threads, "\n")
cat("Anchor workers:", anchor_workers, "\n")
cat("Integrate workers:", integrate_workers, "\n")
cat("Dimensions:", paste(range(dims_use), collapse = "-"), "\n")
cat("Resolutions:", paste(resolution_values, collapse = ", "), "\n")
cat("Save final object:", save_object, "\n\n")

if (!file.exists(bm_file)) stop("Missing BMcounts: ", bm_file)
bm_size <- file.info(bm_file)$size
if (is.na(bm_size) || bm_size != expected_bm_file_size) {
  stop("BMcounts size check failed. Expected ", expected_bm_file_size, " bytes, found ", bm_size)
}

required_packages <- c("Seurat", "Matrix", "ggplot2", "patchwork", "scales")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing required package: ", pkg)
  cat(pkg, "version:", as.character(packageVersion(pkg)), "\n")
}

derive_subject <- function(title) sub("^(UPN[0-9]+)_.*$", "\\1", title)
derive_timepoint <- function(title) {
  if (grepl("_HD[0-9]*", title)) return("HD")
  if (grepl("baseline", title, ignore.case = TRUE)) return("baseline")
  if (grepl("_3M", title)) return("3M")
  if (grepl("_6M", title)) return("6M")
  "unknown"
}
derive_disease <- function(title) ifelse(grepl("_HD[0-9]*", title), "HD", "SAA")
group_label <- function(disease, timepoint) {
  ifelse(disease == "HD", "HD", paste0("SAA_", timepoint))
}

bm_titles <- c(
  "UPN10_BM_GEM_baseline",
  "UPN10_BM_GEM_6M",
  "UPN111_BM_GEM_HD1",
  "UPN112_BM_GEM_HD2",
  "UPN113_BM_GEM_HD3",
  "UPN114_BM_GEM_HD4",
  "UPN11_BM_GEM_baseline",
  "UPN11_BM_6M",
  "UPN12_BM_GEM_baseline",
  "UPN12_BM_6M",
  "UPN13_BM_GEM_baseline",
  "UPN13_BM_6M",
  "UPN14_BM_GEM_baseline",
  "UPN14_BM_3M",
  "UPN14_BM_6M",
  "UPN15_BM_GEM_baseline",
  "UPN15_BM_3M",
  "UPN16_BM_GEM_baseline",
  "UPN16_BM_6M",
  "UPN17_BM_GEM_baseline",
  "UPN17_BM_6M",
  "UPN18_BM_GEM_baseline",
  "UPN18_BM_6M",
  "UPN19_BM_GEM_baseline",
  "UPN19_BM_6M",
  "UPN1_BM_GEM_baseline",
  "UPN1_BM_6M",
  "UPN20_BM_GEM_baseline",
  "UPN20_BM_6M",
  "UPN2_BM_GEM_baseline",
  "UPN2_BM_6M",
  "UPN3_BM_GEM_baseline",
  "UPN3_BM_6M",
  "UPN4_BM_GEM_baseline_repeat",
  "UPN4_BM_GEM_baseline",
  "UPN5_BM_GEM_baseline",
  "UPN6_BM_GEM_baseline",
  "UPN6_BM_3M",
  "UPN7_BM_GEM_baseline",
  "UPN7_BM_6M",
  "UPN8_BM_GEM_baseline",
  "UPN8_BM_6M",
  "UPN9_BM_GEM_baseline",
  "UPN9_BM_3M",
  "UPN9_BM_6M"
)

bm_meta <- data.frame(
  bm_order = seq_along(bm_titles),
  count_suffix = paste0("_", seq_along(bm_titles)),
  geo_accession = paste0("GSM", 7892385 + seq_along(bm_titles)),
  title = bm_titles,
  subject = vapply(bm_titles, derive_subject, character(1)),
  disease = vapply(bm_titles, derive_disease, character(1)),
  timepoint = vapply(bm_titles, derive_timepoint, character(1)),
  is_repeat = grepl("repeat", bm_titles, ignore.case = TRUE),
  stringsAsFactors = FALSE
)
bm_meta$group <- group_label(bm_meta$disease, bm_meta$timepoint)
bm_meta$group <- factor(bm_meta$group, levels = c("HD", "SAA_baseline", "SAA_3M", "SAA_6M"))
write.csv(bm_meta, file.path(out_dir, "97_BM_official_style_sample_metadata.csv"), row.names = FALSE, quote = TRUE)
cat("BM sample group counts:\n")
print(table(bm_meta$group, useNA = "ifany"))

if (dry_run) {
  cat("Dry run complete. Required packages and BMcounts file are available.\n")
  quit(save = "no", status = 0)
}

set_future_plan <- function(workers, label) {
  if (!requireNamespace("future", quietly = TRUE)) return(invisible(FALSE))
  options(future.globals.maxSize = 500 * 1024^3)
  workers <- max(1, as.integer(workers))
  if (workers <= 1) {
    future::plan("sequential")
    cat("Future plan for", label, ": sequential\n")
  } else {
    future::plan("multicore", workers = workers)
    cat("Future plan for", label, ": multicore workers =", workers, "\n")
  }
  invisible(TRUE)
}

set_future_plan(n_threads, "pre-integration preprocessing")

cat("\nLoading BM sparse matrix...\n")
load_env <- new.env(parent = emptyenv())
outer_con <- gzfile(bm_file, "rb")
inner_con <- gzcon(outer_con)
loaded_objects <- load(inner_con, envir = load_env)
close(inner_con)
counts <- load_env[[loaded_objects[1]]]
rm(load_env)
if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
rownames(counts) <- toupper(rownames(counts))
cat("BM matrix class:", class(counts)[1], "\n")
cat("BM matrix dim:", paste(dim(counts), collapse = " x "), "\n")
cat("BM matrix nnzero:", length(counts@x), "\n")
cat("BM matrix object size GB:", round(as.numeric(object.size(counts)) / 1024^3, 3), "\n\n")

cell_ids <- colnames(counts)
count_suffix <- sub("^.*_([0-9]+)$", "_\\1", cell_ids)
count_suffix[!grepl("^_[0-9]+$", count_suffix)] <- NA_character_
if (anyNA(count_suffix)) stop("Failed to parse count suffix from some cell barcodes.")
if (length(setdiff(unique(count_suffix), bm_meta$count_suffix)) > 0) {
  stop("BMcounts contains suffixes not present in metadata.")
}

cell_meta <- bm_meta[match(count_suffix, bm_meta$count_suffix), , drop = FALSE]
rownames(cell_meta) <- cell_ids
cell_meta$orig.ident <- cell_meta$title
cell_meta$group <- factor(cell_meta$group, levels = c("HD", "SAA_baseline", "SAA_3M", "SAA_6M"))

marker_sets <- list(
  HSPC = c("CD34", "KIT", "PROM1", "AVP", "GATA2", "MECOM", "MEIS1", "HLF", "SPINK2"),
  Erythroid = c("HBB", "HBA1", "HBA2", "ALAS2", "GYPA", "KLF1", "GATA1", "AHSP"),
  Megakaryocyte = c("PPBP", "PF4", "ITGA2B", "GP9", "MPL", "VWF"),
  Myeloid = c("LYZ", "LST1", "S100A8", "S100A9", "FCGR3A", "CD14", "MS4A7", "C1QA", "C1QB", "MPO", "ELANE", "AZU1", "CTSG", "FCER1G"),
  T_NK = c("CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", "NKG7", "GNLY", "PRF1", "GZMB", "KLRD1", "IL7R", "CCR7"),
  B_Plasma = c("MS4A1", "CD79A", "CD79B", "CD74", "BANK1", "MZB1", "JCHAIN", "IGHM", "IGKC", "SDC1", "XBP1"),
  Stromal_Endothelial = c("COL1A1", "COL1A2", "DCN", "LUM", "PECAM1", "VWF", "KDR", "CXCL12", "LEPR", "PDGFRA")
)

cat("Assigning broad BM cell types by marker-panel scores...\n")
marker_genes <- unique(unlist(marker_sets, use.names = FALSE))
present_marker_genes <- intersect(marker_genes, rownames(counts))
marker_expr <- counts[present_marker_genes, , drop = FALSE]
lib_size <- Matrix::colSums(counts)
lib_size[!is.finite(lib_size) | lib_size <= 0] <- 1
scale_by_cell <- 10000 / lib_size
marker_expr@x <- log1p(marker_expr@x * rep.int(scale_by_cell, diff(marker_expr@p)))
score_set <- function(genes, mat) {
  genes_present <- intersect(genes, rownames(mat))
  if (length(genes_present) < 1) return(rep(NA_real_, ncol(mat)))
  as.numeric(Matrix::colMeans(mat[genes_present, , drop = FALSE]))
}
marker_scores <- as.matrix(sapply(marker_sets, score_set, mat = marker_expr))
best_idx <- max.col(marker_scores, ties.method = "first")
best_score <- marker_scores[cbind(seq_len(nrow(marker_scores)), best_idx)]
second_score <- apply(marker_scores, 1, function(x) sort(x, decreasing = TRUE)[2])
score_margin <- best_score - second_score
cell_meta$broad_cell_type <- colnames(marker_scores)[best_idx]
cell_meta$broad_cell_type[best_score < 0.02 | score_margin < 0.005] <- "Unknown"
cell_meta$broad_cell_type <- factor(
  cell_meta$broad_cell_type,
  levels = c(names(marker_sets), "Unknown")
)
celltype_summary <- as.data.frame(table(cell_meta$group, cell_meta$broad_cell_type), stringsAsFactors = FALSE)
colnames(celltype_summary) <- c("group", "broad_cell_type", "n_cells")
write.csv(celltype_summary, file.path(out_dir, "98_BM_official_style_celltype_counts_by_group.csv"), row.names = FALSE, quote = TRUE)
rm(marker_expr, marker_scores, scale_by_cell, lib_size)
gc()

cat("\nBuilding per-sample Seurat objects from BMcounts...\n")
seuratObjList_BM <- vector("list", length(bm_titles))
names(seuratObjList_BM) <- bm_titles
for (i in seq_along(bm_titles)) {
  suffix <- paste0("_", i)
  idx <- which(count_suffix == suffix)
  cat("  sample", i, bm_titles[i], "cells:", length(idx), "\n")
  seuratObjList_BM[[i]] <- CreateSeuratObject(
    counts = counts[, idx, drop = FALSE],
    project = bm_titles[i],
    meta.data = cell_meta[idx, , drop = FALSE],
    min.cells = 0,
    min.features = 0
  )
}
rm(counts)
gc()

cat("\nOfficial-style per-sample NormalizeData and FindVariableFeatures...\n")
seuratObjList_BM <- lapply(seuratObjList_BM, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  x
})
features <- SelectIntegrationFeatures(object.list = seuratObjList_BM)
cat("Integration features:", length(features), "\n")

cat("Official-style per-sample ScaleData and RunPCA...\n")
seuratObjList_BM <- lapply(seuratObjList_BM, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, npcs = max(dims_use), verbose = FALSE)
  x
})

reference_idx <- which(bm_meta$disease == "HD")
cat("HD reference indices for RPCA integration:", paste(reference_idx, collapse = ", "), "\n")
cat("Running FindIntegrationAnchors with RPCA...\n")
set_future_plan(anchor_workers, "FindIntegrationAnchors")
sAA.anchors <- FindIntegrationAnchors(
  object.list = seuratObjList_BM,
  reference = reference_idx,
  reduction = "rpca",
  dims = dims_use
)
cat("Running IntegrateData...\n")
set_future_plan(integrate_workers, "IntegrateData")
sAA.combined <- IntegrateData(anchorset = sAA.anchors, dims = dims_use)
DefaultAssay(sAA.combined) <- "integrated"
rm(seuratObjList_BM, sAA.anchors)
gc()

cat("Running official-style ScaleData, PCA, UMAP, neighbors, clusters...\n")
set_future_plan(n_threads, "post-integration PCA/UMAP/clustering")
sAA.combined <- ScaleData(sAA.combined, verbose = FALSE)
sAA.combined <- RunPCA(sAA.combined, npcs = max(dims_use), verbose = FALSE)
sAA.combined <- RunUMAP(sAA.combined, dims = dims_use, verbose = TRUE)
sAA.combined <- FindNeighbors(sAA.combined, dims = dims_use, verbose = TRUE)
for (res in resolution_values) {
  cat("FindClusters resolution:", res, "\n")
  sAA.combined <- FindClusters(sAA.combined, resolution = res, verbose = TRUE)
}

DefaultAssay(sAA.combined) <- "RNA"
meta_out <- cbind(
  as.data.frame(Embeddings(sAA.combined, "umap")),
  sAA.combined@meta.data[, c("count_suffix", "title", "geo_accession", "subject", "disease", "timepoint", "group", "is_repeat", "broad_cell_type"), drop = FALSE]
)
colnames(meta_out)[1:2] <- c("UMAP_1", "UMAP_2")
write.csv(meta_out, file.path(out_dir, "99_BM_official_style_UMAP_cell_metadata.csv"), row.names = TRUE, quote = TRUE)

plot_df <- meta_out
plot_df$group <- factor(plot_df$group, levels = c("HD", "SAA_baseline", "SAA_3M", "SAA_6M"))
plot_df$broad_cell_type <- factor(plot_df$broad_cell_type, levels = c(names(marker_sets), "Unknown"))
celltype_colors <- c(
  HSPC = "#7C3AED",
  Erythroid = "#D7263D",
  Megakaryocyte = "#F59E0B",
  Myeloid = "#2A9D8F",
  T_NK = "#2563EB",
  B_Plasma = "#06B6D4",
  Stromal_Endothelial = "#8D99AE",
  Unknown = "#D1D5DB"
)

theme_umap <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.title = element_blank()
  )

p_all <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = broad_cell_type)) +
  geom_point(size = 0.05, alpha = 0.45, stroke = 0) +
  scale_color_manual(values = celltype_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(title = "GSE247531 BM single cells: broad cell compartments") +
  theme_umap
ggsave(file.path(out_dir, "100_BM_official_style_UMAP_by_broad_celltype.png"), p_all, width = 8.5, height = 6.5, dpi = 300, limitsize = FALSE)
ggsave(file.path(out_dir, "100_BM_official_style_UMAP_by_broad_celltype.pdf"), p_all, width = 8.5, height = 6.5, limitsize = FALSE)

p_split <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = broad_cell_type)) +
  geom_point(size = 0.04, alpha = 0.45, stroke = 0) +
  scale_color_manual(values = celltype_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  facet_wrap(~ group, nrow = 1) +
  labs(title = "GSE247531 BM cell compartments split by group") +
  theme_umap
ggsave(file.path(out_dir, "101_BM_official_style_UMAP_split_by_group_broad_celltype.png"), p_split, width = 16, height = 4.5, dpi = 300, limitsize = FALSE)
ggsave(file.path(out_dir, "101_BM_official_style_UMAP_split_by_group_broad_celltype.pdf"), p_split, width = 16, height = 4.5, limitsize = FALSE)

prop_df <- as.data.frame(table(plot_df$group, plot_df$broad_cell_type), stringsAsFactors = FALSE)
colnames(prop_df) <- c("group", "broad_cell_type", "n_cells")
prop_df <- prop_df[prop_df$group != "" & !is.na(prop_df$group), , drop = FALSE]
prop_df$group_total <- ave(prop_df$n_cells, prop_df$group, FUN = sum)
prop_df$proportion <- prop_df$n_cells / prop_df$group_total
write.csv(prop_df, file.path(out_dir, "102_BM_official_style_celltype_proportion_by_group.csv"), row.names = FALSE, quote = TRUE)
p_prop <- ggplot(prop_df, aes(group, proportion, fill = broad_cell_type)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = celltype_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "BM cell-compartment composition by group", x = NULL, y = "Cell proportion") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.title = element_blank())
ggsave(file.path(out_dir, "103_BM_official_style_celltype_proportion_by_group.png"), p_prop, width = 8.5, height = 5.2, dpi = 300, limitsize = FALSE)
ggsave(file.path(out_dir, "103_BM_official_style_celltype_proportion_by_group.pdf"), p_prop, width = 8.5, height = 5.2, limitsize = FALSE)

if (save_object) {
  cat("Saving final Seurat object RDS...\n")
  saveRDS(sAA.combined, file.path(out_dir, "104_BM_official_style_Seurat_integrated_final.rds"), compress = "xz")
}

cat("\nBM official-style Seurat RPCA UMAP finished at:", as.character(Sys.time()), "\n")
cat("Outputs written to:", out_dir, "\n")
