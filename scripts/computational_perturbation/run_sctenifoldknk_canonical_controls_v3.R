#!/usr/bin/env Rscript

# AA scTenifoldKnk canonical-matrix matched-control extension, version 3.
# This file is additive and does not alter the earlier candidate or control runs.
# All candidate and control runs use one frozen 2328-cell x 3000-gene matrix and
# the same network seed. Before perturbation outcomes are inspected, 20 unique
# nearest controls are frozen for each of the 10 predefined candidates.
# Inferential calibration is performed only by a separate prespecified script
# after matching-quality and support gates have passed.

suppressWarnings(args_raw <- commandArgs(trailingOnly = TRUE))

parse_args <- function(x) {
  out <- list(
    stage = "probe",
    project_root = "PROJECT_ROOT",
    out_root = "",
    candidate_pool = "",
    cell_index = "",
    module_gene_file = "",
    threads = "4",
    common_seed = "20260724",
    shard_id = "1",
    n_shards = "1"
  )
  for (a in x) {
    if (!grepl("^--[^=]+=", a)) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    if (length(kv) == 2L) out[[kv[1]]] <- kv[2]
  }
  out
}

args <- parse_args(args_raw)
stage <- args$stage
out_root <- args$out_root
candidate_pool_path <- args$candidate_pool
cell_index_path <- args$cell_index
module_gene_file <- args$module_gene_file
threads <- max(1L, as.integer(args$threads))
common_seed <- as.integer(args$common_seed)
shard_id <- as.integer(args$shard_id)
n_shards <- as.integer(args$n_shards)

if (!stage %in% c("probe", "design", "benchmark", "candidates",
                  "controls", "full", "summarize")) {
  stop("Unknown --stage: ", stage)
}
if (!nzchar(out_root)) stop("--out_root is required")
if (!nzchar(candidate_pool_path)) stop("--candidate_pool is required")
if (!nzchar(cell_index_path)) stop("--cell_index is required")
if (!is.finite(common_seed)) stop("--common_seed must be an integer")
if (!is.finite(shard_id) || !is.finite(n_shards) ||
    n_shards < 1L || shard_id < 1L || shard_id > n_shards) {
  stop("--shard_id and --n_shards must satisfy 1 <= shard_id <= n_shards")
}

project_root <- args$project_root
bm_file <- file.path(project_root, "01_raw_data", "GSE247531",
                     "GSE247531_BMcounts.Rdata.gz")

core10 <- c("CDK6", "CA2", "PARP1", "KIT", "SYK",
            "GSK3B", "HIF1A", "TOP2A", "TERT", "CD38")
marker_sets <- list(
  HSPC = c("CD34", "KIT", "PROM1", "AVP", "GATA2", "MECOM", "MEIS1",
           "HLF", "SPINK2"),
  Erythroid = c("HBB", "HBA1", "HBA2", "ALAS2", "GYPA", "KLF1",
                "GATA1", "AHSP"),
  Megakaryocyte = c("PPBP", "PF4", "ITGA2B", "GP9", "MPL", "VWF"),
  Myeloid = c("LYZ", "LST1", "S100A8", "S100A9", "FCGR3A", "CD14",
              "MS4A7", "C1QA", "C1QB", "MPO", "ELANE", "AZU1", "CTSG",
              "FCER1G"),
  T_NK = c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1", "CD2",
           "IL7R"),
  B_Plasma = c("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHG1",
               "BANK1"),
  Stromal_Endothelial = c("COL1A1", "COL1A2", "CXCL12", "LEPR",
                          "PECAM1", "VWF", "KDR", "ENG")
)
frozen_markers <- unique(unlist(marker_sets, use.names = FALSE))

required_inputs <- c(bm_file, candidate_pool_path, cell_index_path)
if (nzchar(module_gene_file)) required_inputs <- c(required_inputs, module_gene_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (stage == "probe") {
  checks <- data.frame(
    item = c("BMcounts", "candidate_pool", "cell_index", "module_gene_file",
             "package_Matrix", "package_data.table", "package_scTenifoldNet",
             "package_scTenifoldKnk", "output_root_absent"),
    value = c(
      file.exists(bm_file),
      file.exists(candidate_pool_path),
      file.exists(cell_index_path),
      !nzchar(module_gene_file) || file.exists(module_gene_file),
      requireNamespace("Matrix", quietly = TRUE),
      requireNamespace("data.table", quietly = TRUE),
      requireNamespace("scTenifoldNet", quietly = TRUE),
      requireNamespace("scTenifoldKnk", quietly = TRUE),
      !dir.exists(out_root) && !file.exists(out_root)
    ),
    stringsAsFactors = FALSE
  )
  print(checks, row.names = FALSE)
  if (any(!checks$value)) quit(save = "no", status = 2L)
  quit(save = "no", status = 0L)
}

if (length(missing_inputs)) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "))
}
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})
if (!requireNamespace("scTenifoldNet", quietly = TRUE)) {
  stop("Package scTenifoldNet is required")
}
if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
  stop("Package scTenifoldKnk is required")
}

read_gene_vector <- function(path, preferred = "GeneSymbol") {
  x <- data.table::fread(path, data.table = FALSE, encoding = "UTF-8",
                         showProgress = FALSE)
  if (!ncol(x)) stop("No columns in gene file: ", path)
  col <- if (preferred %in% names(x)) preferred else names(x)[1]
  genes <- trimws(as.character(x[[col]]))
  unique(genes[nzchar(genes) & !is.na(genes)])
}

candidate126 <- read_gene_vector(candidate_pool_path)
if (length(candidate126) != 126L) {
  stop("Expected exactly 126 frozen candidate genes; found ", length(candidate126))
}
if (!all(core10 %in% candidate126)) {
  stop("Core candidate(s) missing from candidate126: ",
       paste(setdiff(core10, candidate126), collapse = ", "))
}

sha256_file <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256"))
  }
  ans <- suppressWarnings(system2("sha256sum", shQuote(path),
                                  stdout = TRUE, stderr = TRUE))
  if (!length(ans) || grepl("not found", ans[1], ignore.case = TRUE)) {
    stop("SHA-256 requires the R package 'digest' or the sha256sum command")
  }
  strsplit(ans[1], "[[:space:]]+")[[1]][1]
}

write_csv_new <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE,
                   fileEncoding = "UTF-8")
}

write_rds_new <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  saveRDS(x, path, compress = "gzip")
}

write_text_new <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  writeLines(x, path, useBytes = TRUE)
}

safe_id <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

load_bm_counts <- function(path) {
  e <- new.env(parent = emptyenv())
  outer <- gzfile(path, "rb")
  inner <- gzcon(outer)
  on.exit(try(close(inner), silent = TRUE), add = TRUE)
  loaded <- load(inner, envir = e)
  if (!"counts" %in% loaded) {
    stop("Object 'counts' was not found in BMcounts; found: ",
         paste(loaded, collapse = ", "))
  }
  counts <- get("counts", envir = e)
  if (!inherits(counts, "Matrix")) counts <- Matrix::Matrix(counts, sparse = TRUE)
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("BMcounts requires gene row names and cell column names")
  }
  counts
}

scale_metric <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

make_wt <- function(dense, seed) {
  set.seed(seed)
  nets <- scTenifoldNet::makeNetworks(
    X = dense, q = 0.9, nNet = 5,
    nCells = min(500L, ncol(dense)),
    scaleScores = TRUE, symmetric = FALSE, nComp = 3, nCores = threads
  )
  td <- scTenifoldNet::tensorDecomposition(
    xList = nets, K = 3, maxError = 1e-5,
    maxIter = 1000, nDecimal = 3
  )
  wt <- as.matrix(scTenifoldKnk:::strictDirection(td$X, lambda = 0))
  diag(wt) <- 0
  wt <- t(wt)
  if (is.null(rownames(wt))) rownames(wt) <- rownames(dense)
  if (is.null(colnames(wt))) colnames(wt) <- rownames(dense)
  wt
}

build_canonical <- function() {
  counts <- load_bm_counts(bm_file)
  idx <- readRDS(cell_index_path)
  if (!is.list(idx) || is.null(idx$selected_cell_ids)) {
    stop("cell_index must contain selected_cell_ids")
  }
  cells <- as.character(idx$selected_cell_ids)
  if (length(cells) != 2328L || anyDuplicated(cells)) {
    stop("Canonical cell index must contain exactly 2328 unique cell IDs; found ",
         length(cells))
  }
  if (!all(cells %in% colnames(counts))) {
    stop("Canonical cell IDs missing from BMcounts: ",
         paste(head(setdiff(cells, colnames(counts)), 10), collapse = ", "))
  }
  if (!all(candidate126 %in% rownames(counts))) {
    stop("All 126 candidates must remain in the canonical network matrix. ",
         "Missing from BMcounts: ",
         paste(setdiff(candidate126, rownames(counts)), collapse = ", "))
  }
  if (!all(core10 %in% rownames(counts))) {
    stop("Core10 missing from BMcounts: ",
         paste(setdiff(core10, rownames(counts)), collapse = ", "))
  }

  sub_all <- counts[, cells, drop = FALSE]
  lib <- Matrix::colSums(sub_all)
  if (any(!is.finite(lib) | lib <= 0)) {
    stop("Nonpositive or nonfinite library size in canonical cells")
  }
  scale_by_cell <- 10000 / lib
  norm <- sub_all
  norm@x <- log1p(norm@x * rep.int(scale_by_cell, diff(norm@p)))
  gene_detect <- Matrix::rowSums(sub_all > 0)
  gene_mean <- Matrix::rowMeans(norm)
  norm_sq <- norm
  norm_sq@x <- norm_sq@x^2
  gene_var <- Matrix::rowMeans(norm_sq) - gene_mean^2
  gene_var[!is.finite(gene_var)] <- 0

  forced <- unique(c(core10, candidate126,
                     intersect(frozen_markers, rownames(counts))))
  if (length(forced) > 3000L) stop("Forced gene set exceeds 3000 genes")
  eligible <- which(gene_detect >= max(5L, ceiling(0.01 * length(cells))))
  ordered <- rownames(counts)[eligible[
    order(-gene_var[eligible], rownames(counts)[eligible])
  ]]
  variable <- setdiff(ordered, forced)
  genes <- c(forced, head(variable, 3000L - length(forced)))
  if (length(genes) != 3000L || anyDuplicated(genes)) {
    stop("Failed to construct exactly 3000 unique canonical genes")
  }
  if (!all(candidate126 %in% genes) || !all(core10 %in% genes)) {
    stop("Candidate or core gene was lost from the canonical panel")
  }

  canonical_sparse <- counts[genes, cells, drop = FALSE]
  canonical_dense <- as.matrix(canonical_sparse)
  storage.mode(canonical_dense) <- "double"
  canonical_norm <- norm[genes, , drop = FALSE]
  metrics <- data.frame(
    gene = genes,
    mean_log1p_cp10k = as.numeric(Matrix::rowMeans(canonical_norm)),
    detected_fraction = as.numeric(Matrix::rowMeans(canonical_sparse > 0)),
    variance_log1p_cp10k = as.numeric(gene_var[genes]),
    in_core10 = genes %in% core10,
    in_candidate126 = genes %in% candidate126,
    in_frozen_marker = genes %in% frozen_markers,
    stringsAsFactors = FALSE
  )

  wt <- make_wt(canonical_dense, common_seed)
  if (!all(dim(wt) == c(3000L, 3000L))) {
    stop("Unexpected baseline network dimensions: ",
         paste(dim(wt), collapse = "x"))
  }
  metrics$out_degree <- rowSums(abs(wt) > 0)
  metrics$out_strength <- rowSums(abs(wt))
  metrics$z_mean <- scale_metric(metrics$mean_log1p_cp10k)
  frac <- pmin(1 - 1e-6, pmax(1e-6, metrics$detected_fraction))
  metrics$z_detect <- scale_metric(qlogis(frac))
  metrics$z_degree <- scale_metric(log1p(metrics$out_degree))
  metrics$z_strength <- scale_metric(log1p(metrics$out_strength))
  metrics$eligible_control <- (
    !metrics$in_candidate126 &
      metrics$detected_fraction >= 0.02 &
      metrics$mean_log1p_cp10k > 0 &
      metrics$out_degree > 0 &
      metrics$out_strength > 0
  )

  list(
    counts = canonical_dense,
    cells = cells,
    genes = genes,
    metrics = metrics,
    forced_genes = forced,
    common_seed = common_seed,
    source = list(
      bm_file = bm_file,
      bm_sha256 = sha256_file(bm_file),
      candidate_pool = candidate_pool_path,
      candidate_pool_sha256 = sha256_file(candidate_pool_path),
      cell_index = cell_index_path,
      cell_index_sha256 = sha256_file(cell_index_path)
    )
  )
}

make_design <- function(canonical) {
  metrics <- canonical$metrics
  used <- character()
  matched_rows <- list()
  for (target in core10) {
    trow <- metrics[metrics$gene == target, , drop = FALSE]
    if (nrow(trow) != 1L) stop("Target missing or duplicated: ", target)
    pool <- metrics[metrics$eligible_control & !metrics$gene %in% used,
                    , drop = FALSE]
    pool$match_distance <- sqrt(
      (pool$z_mean - trow$z_mean)^2 +
        (pool$z_detect - trow$z_detect)^2 +
        (pool$z_degree - trow$z_degree)^2 +
        (pool$z_strength - trow$z_strength)^2
    )
    pool <- pool[order(pool$match_distance, pool$gene), , drop = FALSE]
    if (nrow(pool) < 20L) stop("Fewer than 20 matched controls for ", target)
    pick <- head(pool, 20L)
    used <- c(used, pick$gene)
    matched_rows[[length(matched_rows) + 1L]] <- data.frame(
      run_role = "matched_control",
      candidate = target,
      gKO = pick$gene,
      control_type = "matched",
      control_rank = seq_len(20L),
      common_seed = common_seed,
      match_distance = pick$match_distance,
      candidate_mean_log1p_cp10k = trow$mean_log1p_cp10k,
      control_mean_log1p_cp10k = pick$mean_log1p_cp10k,
      candidate_detected_fraction = trow$detected_fraction,
      control_detected_fraction = pick$detected_fraction,
      candidate_out_degree = trow$out_degree,
      control_out_degree = pick$out_degree,
      candidate_out_strength = trow$out_strength,
      control_out_strength = pick$out_strength,
      stringsAsFactors = FALSE
    )
  }
  matched <- do.call(rbind, matched_rows)

  candidate <- data.frame(
    run_role = "candidate",
    candidate = core10,
    gKO = core10,
    control_type = "candidate",
    control_rank = NA_integer_,
    common_seed = common_seed,
    match_distance = NA_real_,
    candidate_mean_log1p_cp10k =
      metrics$mean_log1p_cp10k[match(core10, metrics$gene)],
    control_mean_log1p_cp10k = NA_real_,
    candidate_detected_fraction =
      metrics$detected_fraction[match(core10, metrics$gene)],
    control_detected_fraction = NA_real_,
    candidate_out_degree = metrics$out_degree[match(core10, metrics$gene)],
    control_out_degree = NA_real_,
    candidate_out_strength = metrics$out_strength[match(core10, metrics$gene)],
    control_out_strength = NA_real_,
    stringsAsFactors = FALSE
  )
  manifest <- rbind(candidate, matched)
  manifest$run_id <- sprintf(
    "%03d_%s_%s",
    seq_len(nrow(manifest)), manifest$run_role, safe_id(manifest$gKO)
  )
  manifest <- manifest[, c("run_id", setdiff(names(manifest), "run_id"))]
  if (nrow(manifest) != 210L) stop("Expected 210 total runs")
  if (sum(manifest$run_role == "matched_control") != 200L) {
    stop("Expected 200 matched-control runs")
  }
  if (length(unique(manifest$common_seed)) != 1L) {
    stop("All runs must use one common frozen network seed")
  }
  manifest
}

rank_result <- function(x) {
  df <- as.data.frame(x)
  if (!"gene" %in% names(df)) df$gene <- rownames(df)
  for (nm in c("distance", "Z", "FC", "p.value", "p.adj")) {
    if (!nm %in% names(df)) df[[nm]] <- NA_real_
    df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
  }
  df
}

module_genes <- if (nzchar(module_gene_file)) {
  read_gene_vector(module_gene_file)
} else {
  character()
}

endpoint_row <- function(df, manifest_row, runtime_minutes) {
  dr <- rank_result(df)
  dr_eval <- dr[dr$gene != manifest_row$gKO, , drop = FALSE]
  max_z <- if (any(is.finite(dr_eval$Z))) {
    max(abs(dr_eval$Z), na.rm = TRUE)
  } else {
    NA_real_
  }
  module_score <- if (length(module_genes)) {
    sum(abs(dr_eval$Z[dr_eval$gene %in% module_genes]), na.rm = TRUE)
  } else {
    NA_real_
  }
  data.frame(
    run_id = manifest_row$run_id,
    run_role = manifest_row$run_role,
    candidate = manifest_row$candidate,
    gKO = manifest_row$gKO,
    common_seed = manifest_row$common_seed,
    n_sig_excluding_gKO_padj_0_05 =
      sum(dr_eval$p.adj < 0.05, na.rm = TRUE),
    max_abs_Z_excluding_gKO_exploratory = max_z,
    module_abs_Z_sum_excluding_gKO_exploratory = module_score,
    runtime_minutes = runtime_minutes,
    stringsAsFactors = FALSE
  )
}

check_freeze <- function() {
  manifest_path <- file.path(out_root, "design", "run_manifest_210.csv")
  freeze_path <- file.path(out_root, "design", "DESIGN_FREEZE_V3.txt")
  input_path <- file.path(out_root, "inputs", "canonical_2328x3000.rds")
  for (p in c(manifest_path, freeze_path, input_path)) {
    if (!file.exists(p)) stop("Frozen design input missing: ", p)
  }
  freeze <- readLines(freeze_path, warn = FALSE)
  expected <- sub("^manifest_sha256=", "",
                  freeze[grepl("^manifest_sha256=", freeze)])
  if (length(expected) != 1L || sha256_file(manifest_path) != expected) {
    stop("Frozen manifest SHA-256 mismatch")
  }
  list(
    manifest = utils::read.csv(manifest_path, stringsAsFactors = FALSE,
                               check.names = FALSE),
    canonical = readRDS(input_path)
  )
}

run_one <- function(canonical, row) {
  set.seed(as.integer(row$common_seed))
  started <- Sys.time()
  result <- scTenifoldKnk::scTenifoldKnk(
    countMatrix = canonical$counts,
    gKO = row$gKO,
    qc = FALSE,
    nc_nNet = 5,
    nc_nCells = min(500L, ncol(canonical$counts)),
    td_K = 3,
    nCores = threads
  )
  runtime <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  dr <- as.data.frame(result$diffRegulation)
  if (!"gene" %in% names(dr)) dr$gene <- rownames(dr)
  run_dir <- file.path(out_root, "runs", row$run_role)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(run_dir, paste0(row$run_id, "_diffRegulation.csv"))
  endpoint_path <- file.path(run_dir, paste0(row$run_id, "_endpoint.csv"))
  if (file.exists(result_path) || file.exists(endpoint_path)) {
    stop("Refusing to overwrite existing run output: ", row$run_id)
  }
  write_csv_new(dr, result_path)
  ep <- endpoint_row(dr, row, runtime)
  write_csv_new(ep, endpoint_path)
  ep
}

if (stage == "design") {
  if (dir.exists(out_root) || file.exists(out_root)) {
    stop("Refusing to overwrite existing output root: ", out_root)
  }
  dir.create(file.path(out_root, "design"), recursive = TRUE)
  dir.create(file.path(out_root, "inputs"), recursive = TRUE)
  dir.create(file.path(out_root, "runs"), recursive = TRUE)
  dir.create(file.path(out_root, "summary"), recursive = TRUE)
  canonical <- build_canonical()
  manifest <- make_design(canonical)
  canonical_path <- file.path(out_root, "inputs", "canonical_2328x3000.rds")
  metrics_path <- file.path(out_root, "design", "canonical_gene_metrics.csv")
  manifest_path <- file.path(out_root, "design", "run_manifest_210.csv")
  write_rds_new(canonical, canonical_path)
  write_csv_new(canonical$metrics, metrics_path)
  write_csv_new(manifest, manifest_path)
  write_text_new(
    c(
      "version=3",
      "matrix=2328_cells_x_3000_genes",
      paste0("common_seed=", common_seed),
      "core_candidates=10",
      "candidate_pool_retained_in_matrix=126",
      "candidate_pool_excluded_from_controls=126",
      "matched_controls_per_candidate=20",
      "matched_control_runs=200",
      "shared_random_control_runs=0",
      "candidate_runs=10",
      "total_runs=210",
      paste0("canonical_rds_sha256=", sha256_file(canonical_path)),
      paste0("metrics_sha256=", sha256_file(metrics_path)),
      paste0("manifest_sha256=", sha256_file(manifest_path)),
      "inference=deferred_to_prespecified_covariate_adjusted_calibration",
      "empirical_P=only_after_matching_and_support_gates",
      "BH=only_for_candidates_passing_prespecified_gates"
    ),
    file.path(out_root, "design", "DESIGN_FREEZE_V3.txt")
  )
  cat("DESIGN_V3_COMPLETE\n", manifest_path, "\n")
  quit(save = "no", status = 0L)
}

frozen <- check_freeze()
manifest <- frozen$manifest
canonical <- frozen$canonical
if (!all(dim(canonical$counts) == c(3000L, 2328L))) {
  stop("Frozen canonical matrix is not 3000 genes x 2328 cells")
}
if (!all(candidate126 %in% rownames(canonical$counts))) {
  stop("Frozen canonical matrix does not retain all 126 candidates")
}
if (length(unique(manifest$common_seed)) != 1L ||
    unique(manifest$common_seed) != common_seed) {
  stop("Common seed mismatch between invocation and frozen manifest")
}

if (stage == "benchmark") {
  selected <- c(
    which(manifest$run_role == "candidate")[1],
    which(manifest$run_role == "matched_control")[1],
    which(manifest$run_role == "matched_control" &
            manifest$candidate == core10[2])[1]
  )
} else if (stage == "candidates") {
  selected <- which(manifest$run_role == "candidate")
} else if (stage == "controls") {
  selected <- which(manifest$run_role != "candidate")
} else if (stage == "full") {
  selected <- seq_len(nrow(manifest))
} else if (stage == "summarize") {
  endpoint_files <- list.files(file.path(out_root, "runs"),
                               pattern = "_endpoint[.]csv$",
                               recursive = TRUE, full.names = TRUE)
  if (!length(endpoint_files)) stop("No endpoint files to summarize")
  endpoints <- do.call(rbind, lapply(endpoint_files, function(p) {
    utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  }))
  endpoints <- endpoints[order(match(endpoints$run_id, manifest$run_id)), ]
  write_csv_new(
    endpoints,
    file.path(out_root, "summary", "descriptive_endpoints_all_runs.csv")
  )
  cat("DESCRIPTIVE_SUMMARY_COMPLETE\n")
  quit(save = "no", status = 0L)
}

if (stage != "summarize" && n_shards > 1L) {
  selected <- selected[((seq_along(selected) - 1L) %% n_shards) ==
                         (shard_id - 1L)]
}

for (i in selected) {
  row <- manifest[i, , drop = FALSE]
  result_path <- file.path(out_root, "runs", row$run_role,
                           paste0(row$run_id, "_diffRegulation.csv"))
  endpoint_path <- file.path(out_root, "runs", row$run_role,
                             paste0(row$run_id, "_endpoint.csv"))
  if (file.exists(result_path) && file.exists(endpoint_path)) {
    cat("SKIP_COMPLETE", row$run_id, "\n")
    next
  }
  if (file.exists(result_path) || file.exists(endpoint_path)) {
    stop("Partial pre-existing run requires manual review: ", row$run_id)
  }
  ep <- run_one(canonical, row)
  cat("RUN_COMPLETE", row$run_id,
      sprintf("runtime_minutes=%.2f", ep$runtime_minutes), "\n")
  gc()
}
cat(toupper(stage), "_V3_COMPLETE shard=", shard_id, "/", n_shards, "\n",
    sep = "")
