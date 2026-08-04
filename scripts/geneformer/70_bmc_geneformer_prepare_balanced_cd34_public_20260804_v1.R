#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(paste(
    "Usage: Rscript 70_bmc_geneformer_prepare_balanced_cd34_public_20260804_v1.R",
    "GSE247531_CD34counts.Rdata.gz DESIGN_TABLE.tsv OUTPUT_DIR"
  ))
}

input_counts <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
design_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[3]], winslash = "/", mustWork = FALSE)
seed <- 20260802L
cells_per_subject_timepoint <- 64L

if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Refusing to write into a non-empty output directory: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# The archived GEO file contains two gzip layers. Open each layer explicitly;
# no decompressed copy is written and the source archive remains read-only.
outer <- gzfile(input_counts, open = "rb")
inner <- gzcon(outer)
env <- new.env(parent = emptyenv())
loaded <- load(inner, envir = env)
close(inner)
if (!identical(loaded, "counts")) {
  stop("Unexpected RData objects: ", paste(loaded, collapse = ","))
}
counts <- get("counts", envir = env)
if (!inherits(counts, "dgCMatrix")) {
  stop("Expected dgCMatrix; observed: ", paste(class(counts), collapse = "/"))
}

design <- read.delim(design_file, check.names = FALSE, stringsAsFactors = FALSE)
required_design <- c("count_suffix", "geo_accession", "subject", "disease", "timepoint", "is_repeat")
if (!all(required_design %in% names(design))) {
  stop("Design table lacks required fields: ", paste(setdiff(required_design, names(design)), collapse = ","))
}
design$count_suffix <- as.character(design$count_suffix)
if (anyDuplicated(design$count_suffix)) stop("Duplicate count_suffix values in design table")

cell_id <- colnames(counts)
count_suffix <- sub("^.*(_[0-9]+)$", "\\1", cell_id)
design_index <- match(count_suffix, design$count_suffix)
if (anyNA(design_index)) {
  stop("Unmapped cell suffixes: ", paste(head(unique(count_suffix[is.na(design_index)]), 20), collapse = ","))
}
cell_meta <- design[design_index, required_design, drop = FALSE]
cell_meta$cell_id <- cell_id
cell_meta$count_suffix <- count_suffix
cell_meta$timepoint[is.na(cell_meta$timepoint) | cell_meta$timepoint == ""] <- "HD"
cell_meta$subject_timepoint <- paste(cell_meta$subject, cell_meta$timepoint, sep = "__")

set.seed(seed)
split_indices <- split(seq_len(ncol(counts)), cell_meta$subject_timepoint)
sampled_indices <- unlist(lapply(split_indices, function(idx) {
  n_take <- min(length(idx), cells_per_subject_timepoint)
  sort(sample(idx, size = n_take, replace = FALSE))
}), use.names = FALSE)
sampled_indices <- sort(sampled_indices)

sampled <- counts[, sampled_indices, drop = FALSE]
sampled_meta <- cell_meta[sampled_indices, , drop = FALSE]
sampled_meta$sampled_order <- seq_len(nrow(sampled_meta))
sampled_meta$n_counts <- as.numeric(colSums(sampled))
sampled_meta$filter_pass <- 1L
if (any(sampled_meta$n_counts <= 0)) stop("Sampled cells with zero total counts detected")

matrix_file <- file.path(output_dir, "GSE247531_CD34_balanced_cells_by_genes_mvp_v1.mtx")
genes_file <- file.path(output_dir, "GSE247531_CD34_balanced_gene_symbols_mvp_v1.tsv")
cells_file <- file.path(output_dir, "GSE247531_CD34_balanced_cell_metadata_mvp_v1.tsv")
group_file <- file.path(output_dir, "GSE247531_CD34_balanced_group_counts_mvp_v1.tsv")
manifest_file <- file.path(output_dir, "GSE247531_CD34_balanced_input_manifest_mvp_v1.tsv")

writeMM(t(sampled), matrix_file)
write.table(data.frame(gene_symbol = rownames(sampled), stringsAsFactors = FALSE), genes_file,
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sampled_meta, cells_file, sep = "\t", quote = FALSE, row.names = FALSE)

group_counts <- as.data.frame(table(subject = sampled_meta$subject,
                                    disease = sampled_meta$disease,
                                    timepoint = sampled_meta$timepoint), stringsAsFactors = FALSE)
group_counts <- group_counts[group_counts$Freq > 0, , drop = FALSE]
write.table(group_counts, group_file, sep = "\t", quote = FALSE, row.names = FALSE)

manifest <- data.frame(
  key = c("source_counts_basename", "source_design_basename", "archive_layers", "seed",
          "cells_per_subject_timepoint", "source_dimensions", "sampled_dimensions",
          "sampled_nonzero", "n_subjects", "n_subject_timepoints", "n_HD_subjects", "n_SAA_subjects"),
  value = c(basename(input_counts), basename(design_file), "gzip_x2", seed,
            cells_per_subject_timepoint, paste(dim(counts), collapse = "x"),
            paste(dim(sampled), collapse = "x"), length(sampled@x),
            length(unique(sampled_meta$subject)), length(unique(sampled_meta$subject_timepoint)),
            length(unique(sampled_meta$subject[sampled_meta$disease == "HD"])),
            length(unique(sampled_meta$subject[sampled_meta$disease == "SAA"]))),
  stringsAsFactors = FALSE
)
write.table(manifest, manifest_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("PREPARE_OK\n")
cat("SOURCE_DIM\t", paste(dim(counts), collapse = "x"), "\n", sep = "")
cat("SAMPLED_DIM\t", paste(dim(sampled), collapse = "x"), "\n", sep = "")
cat("SAMPLED_NONZERO\t", length(sampled@x), "\n", sep = "")
cat("SUBJECTS\t", length(unique(sampled_meta$subject)), "\n", sep = "")
cat("SUBJECT_TIMEPOINTS\t", length(unique(sampled_meta$subject_timepoint)), "\n", sep = "")
cat("OUTPUT_DIR\t", output_dir, "\n", sep = "")
