#!/usr/bin/env Rscript

# Prespecified calibration for the additive 210-run scTenifoldKnk v3 batch.
# This script never changes the frozen candidate order. Empirical P and BH
# values are emitted only when the declared matching and null-homogeneity gates
# pass; otherwise the corresponding output is reported without calibration.

args_raw <- commandArgs(trailingOnly = TRUE)
parse_args <- function(x) {
  out <- list(out_root = "")
  for (a in x) {
    if (!grepl("^--[^=]+=", a)) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    if (length(kv) == 2L) out[[kv[1]]] <- kv[2]
  }
  out
}
args <- parse_args(args_raw)
out_root <- args$out_root
if (!nzchar(out_root)) stop("--out_root is required")

manifest_path <- file.path(out_root, "design", "run_manifest_210.csv")
endpoint_path <- file.path(
  out_root, "summary", "descriptive_endpoints_all_runs.csv"
)
for (p in c(manifest_path, endpoint_path)) {
  if (!file.exists(p)) stop("Required input missing: ", p)
}

manifest <- read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
endpoints <- read.csv(
  endpoint_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(manifest) != 210L || nrow(endpoints) != 210L) {
  stop("Expected exactly 210 manifest and endpoint rows")
}
if (anyDuplicated(manifest$run_id) || anyDuplicated(endpoints$run_id)) {
  stop("Duplicated run_id detected")
}
if (!setequal(manifest$run_id, endpoints$run_id)) {
  stop("Manifest and endpoint run IDs differ")
}
dat <- merge(
  manifest, endpoints,
  by = c("run_id", "run_role", "candidate", "gKO", "common_seed"),
  all = FALSE, sort = FALSE
)
dat <- dat[match(manifest$run_id, dat$run_id), , drop = FALSE]
if (nrow(dat) != 210L || any(is.na(dat$run_id))) {
  stop("Failed to align all 210 runs")
}

core10 <- c(
  "CDK6", "CA2", "PARP1", "KIT", "SYK",
  "GSK3B", "HIF1A", "TOP2A", "TERT", "CD38"
)
endpoint_col <- "n_sig_excluding_gKO_padj_0_05"
if (!endpoint_col %in% names(dat)) stop("Primary endpoint column missing")

metric_pairs <- list(
  mean_log1p_cp10k = c(
    "candidate_mean_log1p_cp10k", "control_mean_log1p_cp10k"
  ),
  detected_fraction = c(
    "candidate_detected_fraction", "control_detected_fraction"
  ),
  out_degree = c("candidate_out_degree", "control_out_degree"),
  out_strength = c("candidate_out_strength", "control_out_strength")
)

candidate_rows <- list()
null_rows <- list()
for (target in core10) {
  cand <- dat[dat$run_role == "candidate" & dat$gKO == target, , drop = FALSE]
  ctrl <- dat[
    dat$run_role == "matched_control" & dat$candidate == target,
    , drop = FALSE
  ]
  if (nrow(cand) != 1L || nrow(ctrl) != 20L) {
    stop("Expected one candidate and 20 matched controls for ", target)
  }
  y <- as.numeric(ctrl[[endpoint_col]])
  yc <- as.numeric(cand[[endpoint_col]])
  if (any(!is.finite(c(y, yc)))) {
    stop("Nonfinite primary endpoint for ", target)
  }

  balance <- numeric(length(metric_pairs))
  covered <- logical(length(metric_pairs))
  names(balance) <- names(metric_pairs)
  names(covered) <- names(metric_pairs)
  for (nm in names(metric_pairs)) {
    pair <- metric_pairs[[nm]]
    candidate_value <- unique(as.numeric(ctrl[[pair[1]]]))
    control_values <- as.numeric(ctrl[[pair[2]]])
    if (length(candidate_value) != 1L ||
        any(!is.finite(c(candidate_value, control_values)))) {
      balance[nm] <- Inf
      covered[nm] <- FALSE
      next
    }
    s <- sd(control_values)
    balance[nm] <- if (is.finite(s) && s > 0) {
      abs(candidate_value - mean(control_values)) / s
    } else {
      Inf
    }
    covered[nm] <- (
      candidate_value >= min(control_values) &&
        candidate_value <= max(control_values)
    )
  }

  sy <- sd(y)
  candidate_z <- if (is.finite(sy) && sy > 0) {
    (yc - mean(y)) / sy
  } else {
    NA_real_
  }
  local_percentile <- (
    1 + sum(y <= yc, na.rm = TRUE)
  ) / (length(y) + 1)

  loo_z <- rep(NA_real_, length(y))
  for (i in seq_along(y)) {
    ref <- y[-i]
    sref <- sd(ref)
    if (is.finite(sref) && sref > 0) {
      loo_z[i] <- (y[i] - mean(ref)) / sref
    }
  }
  null_rows[[target]] <- data.frame(
    candidate = target,
    control_gene = ctrl$gKO,
    endpoint = y,
    leave_one_out_standardized_residual = loo_z,
    stringsAsFactors = FALSE
  )

  individual_gate <- (
    all(is.finite(balance)) &&
      max(balance) <= 1.0 &&
      all(covered) &&
      is.finite(candidate_z) &&
      sum(is.finite(loo_z)) == 20L
  )
  reasons <- character()
  if (max(balance) > 1.0) reasons <- c(reasons, "matching_imbalance_gt_1SD")
  if (!all(covered)) reasons <- c(reasons, "candidate_outside_control_range")
  if (!is.finite(candidate_z)) reasons <- c(reasons, "zero_endpoint_variance")
  if (sum(is.finite(loo_z)) != 20L) {
    reasons <- c(reasons, "nonfinite_leave_one_out_residual")
  }

  candidate_rows[[target]] <- data.frame(
    candidate = target,
    candidate_endpoint = yc,
    matched_control_mean = mean(y),
    matched_control_median = median(y),
    matched_control_sd = sy,
    endpoint_difference_from_median = yc - median(y),
    candidate_standardized_residual = candidate_z,
    matched_control_local_percentile = local_percentile,
    max_abs_standardized_matching_imbalance = max(balance),
    all_four_metrics_range_covered = all(covered),
    individual_gate_pass = individual_gate,
    individual_gate_reason = if (length(reasons)) {
      paste(reasons, collapse = ";")
    } else {
      "PASS"
    },
    stringsAsFactors = FALSE
  )
}

candidates <- do.call(rbind, candidate_rows)
null <- do.call(rbind, null_rows)
finite_null <- null[is.finite(null$leave_one_out_standardized_residual), ]
if (nrow(finite_null) != 200L) {
  stop("Expected 200 finite leave-one-out null residuals")
}
homogeneity <- kruskal.test(
  leave_one_out_standardized_residual ~ candidate,
  data = finite_null
)
global_gate <- is.finite(homogeneity$p.value) && homogeneity$p.value >= 0.01

candidates$pooled_null_homogeneity_p <- homogeneity$p.value
candidates$global_homogeneity_gate_pass <- global_gate
candidates$empirical_calibration_p <- NA_real_
for (i in seq_len(nrow(candidates))) {
  if (candidates$individual_gate_pass[i] && global_gate) {
    zc <- candidates$candidate_standardized_residual[i]
    candidates$empirical_calibration_p[i] <- (
      1 + sum(
        finite_null$leave_one_out_standardized_residual >= zc
      )
    ) / (nrow(finite_null) + 1)
  }
}
candidates$BH_q <- NA_real_
eligible <- which(is.finite(candidates$empirical_calibration_p))
if (length(eligible)) {
  candidates$BH_q[eligible] <- p.adjust(
    candidates$empirical_calibration_p[eligible], method = "BH"
  )
}
candidates$reporting_status <- ifelse(
  is.finite(candidates$empirical_calibration_p),
  "calibrated_empirical_result",
  "not_calibration_eligible"
)
candidates <- candidates[match(core10, candidates$candidate), , drop = FALSE]

summary_dir <- file.path(out_root, "summary")
candidate_out <- file.path(
  summary_dir, "candidate_matched_calibration_v3.csv"
)
null_out <- file.path(summary_dir, "pooled_matched_null_residuals_v3.csv")
gate_out <- file.path(summary_dir, "calibration_gates_v3.txt")
for (p in c(candidate_out, null_out, gate_out)) {
  if (file.exists(p)) stop("Refusing to overwrite: ", p)
}
write.csv(
  candidates, candidate_out, row.names = FALSE, quote = TRUE,
  fileEncoding = "UTF-8"
)
write.csv(
  null, null_out, row.names = FALSE, quote = TRUE,
  fileEncoding = "UTF-8"
)
writeLines(
  c(
    "version=3",
    "primary_endpoint=n_sig_excluding_gKO_padj_0_05",
    "matched_controls_per_candidate=20",
    "pooled_leave_one_out_null_residuals=200",
    paste0("kruskal_wallis_p=", format(homogeneity$p.value, digits = 12)),
    paste0("global_homogeneity_gate_pass=", global_gate),
    paste0(
      "candidates_passing_individual_gate=",
      sum(candidates$individual_gate_pass)
    ),
    paste0(
      "candidates_with_calibrated_empirical_p=",
      sum(is.finite(candidates$empirical_calibration_p))
    ),
    "claim_boundary=no_causal_or_target_specificity_claim"
  ),
  gate_out, useBytes = TRUE
)
cat("MATCHED_CALIBRATION_V3_COMPLETE\n")
cat(candidate_out, "\n", null_out, "\n", gate_out, "\n")
