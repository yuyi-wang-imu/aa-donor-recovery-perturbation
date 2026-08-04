#!/usr/bin/env Rscript

# Balance-triggered sensitivity calibration for the completed v3 batch.
#
# The global rematching assignment was frozen using baseline covariates only.
# This script applies that locked, disjoint assignment to the already completed
# perturbation endpoints. It does not alter or replace the original v3 primary
# calibration.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: Rscript summarize_sctenifoldknk_global_unique_sensitivity_v1.R ",
    "<manifest.csv> <endpoints.csv> <assignments.csv> <output_directory>"
  )
}
manifest_path <- normalizePath(args[[1]], mustWork = TRUE)
endpoint_path <- normalizePath(args[[2]], mustWork = TRUE)
assignment_path <- normalizePath(args[[3]], mustWork = TRUE)
output_dir <- normalizePath(args[[4]], mustWork = TRUE)

candidate_out <- file.path(
  output_dir, "candidate_global_unique_calibration_sensitivity_v1.csv"
)
null_out <- file.path(
  output_dir, "pooled_global_unique_null_residuals_sensitivity_v1.csv"
)
gate_out <- file.path(
  output_dir, "global_unique_calibration_sensitivity_v1.txt"
)
for (p in c(candidate_out, null_out, gate_out)) {
  if (file.exists(p)) stop("Refusing to overwrite: ", p)
}

manifest <- read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
endpoints <- read.csv(
  endpoint_path, stringsAsFactors = FALSE, check.names = FALSE
)
assignments <- read.csv(
  assignment_path, stringsAsFactors = FALSE, check.names = FALSE
)
core10 <- c(
  "CDK6", "CA2", "PARP1", "KIT", "SYK",
  "GSK3B", "HIF1A", "TOP2A", "TERT", "CD38"
)
endpoint_col <- "n_sig_excluding_gKO_padj_0_05"
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

if (nrow(assignments) != 200L || anyDuplicated(assignments$control_gene)) {
  stop("Assignments must contain 200 unique controls")
}
if (!endpoint_col %in% names(endpoints)) stop("Primary endpoint missing")

dat <- merge(
  manifest, endpoints[, c("run_id", endpoint_col)],
  by = "run_id", all.x = TRUE, sort = FALSE
)
if (nrow(dat) != 210L || any(!is.finite(dat[[endpoint_col]]))) {
  stop("Failed to align all 210 endpoints")
}

candidate_rows <- list()
null_rows <- list()
for (target in core10) {
  cand <- dat[dat$run_role == "candidate" & dat$gKO == target, ,
              drop = FALSE]
  assigned_genes <- assignments$control_gene[
    assignments$candidate == target
  ]
  ctrl <- dat[
    dat$run_role == "matched_control" & dat$gKO %in% assigned_genes, ,
    drop = FALSE
  ]
  if (nrow(cand) != 1L || nrow(ctrl) != 20L) {
    stop("Candidate/control mismatch: ", target)
  }

  y <- as.numeric(ctrl[[endpoint_col]])
  yc <- as.numeric(cand[[endpoint_col]])
  balance <- numeric(length(metric_pairs))
  covered <- logical(length(metric_pairs))
  names(balance) <- names(metric_pairs)
  names(covered) <- names(metric_pairs)
  for (nm in names(metric_pairs)) {
    pair <- metric_pairs[[nm]]
    candidate_value <- as.numeric(cand[[pair[[1L]]]])
    control_values <- as.numeric(ctrl[[pair[[2L]]]])
    s <- sd(control_values)
    balance[[nm]] <- if (is.finite(s) && s > 0) {
      abs(candidate_value - mean(control_values)) / s
    } else {
      Inf
    }
    covered[[nm]] <- candidate_value >= min(control_values) &&
      candidate_value <= max(control_values)
  }

  sy <- sd(y)
  candidate_z <- if (is.finite(sy) && sy > 0) {
    (yc - mean(y)) / sy
  } else {
    NA_real_
  }
  loo_z <- rep(NA_real_, length(y))
  for (i in seq_along(y)) {
    ref <- y[-i]
    sref <- sd(ref)
    if (is.finite(sref) && sref > 0) {
      loo_z[[i]] <- (y[[i]] - mean(ref)) / sref
    }
  }
  strict_gate <- all(is.finite(balance)) &&
    max(balance) <= 1 &&
    all(covered) &&
    is.finite(candidate_z) &&
    sum(is.finite(loo_z)) == 20L
  reasons <- character()
  if (max(balance) > 1) {
    reasons <- c(reasons, "matching_imbalance_gt_1SD")
  }
  if (!all(covered)) {
    reasons <- c(reasons, "candidate_outside_control_range")
  }
  if (!is.finite(candidate_z)) {
    reasons <- c(reasons, "zero_endpoint_variance")
  }
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
    matched_control_local_percentile = (
      1 + sum(y <= yc)
    ) / (length(y) + 1),
    max_abs_standardized_matching_imbalance = max(balance),
    all_four_metrics_range_covered = all(covered),
    individual_gate_pass = strict_gate,
    individual_gate_reason = if (length(reasons)) {
      paste(reasons, collapse = ";")
    } else {
      "PASS"
    },
    stringsAsFactors = FALSE
  )
  null_rows[[target]] <- data.frame(
    candidate = target,
    control_gene = ctrl$gKO,
    endpoint = y,
    leave_one_out_standardized_residual = loo_z,
    eligible_group = strict_gate,
    stringsAsFactors = FALSE
  )
}

candidates <- do.call(rbind, candidate_rows)
null_all <- do.call(rbind, null_rows)
rownames(candidates) <- NULL
rownames(null_all) <- NULL
finite_null <- null_all[
  null_all$eligible_group &
    is.finite(null_all$leave_one_out_standardized_residual),
  , drop = FALSE
]

eligible_groups <- unique(finite_null$candidate)
if (length(eligible_groups) < 2L) {
  stop("Too few eligible groups for null homogeneity test")
}
homogeneity <- kruskal.test(
  leave_one_out_standardized_residual ~ candidate,
  data = finite_null
)
global_gate <- is.finite(homogeneity$p.value) &&
  homogeneity$p.value >= 0.01

candidates$pooled_null_homogeneity_p <- homogeneity$p.value
candidates$global_homogeneity_gate_pass <- global_gate
candidates$empirical_calibration_p <- NA_real_
for (i in seq_len(nrow(candidates))) {
  if (candidates$individual_gate_pass[[i]] && global_gate) {
    zc <- candidates$candidate_standardized_residual[[i]]
    candidates$empirical_calibration_p[[i]] <- (
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
  "calibrated_sensitivity_result",
  "not_calibration_eligible"
)
candidates <- candidates[match(core10, candidates$candidate), , drop = FALSE]

write.csv(candidates, candidate_out, row.names = FALSE, na = "")
write.csv(finite_null, null_out, row.names = FALSE, na = "")
writeLines(
  c(
    "analysis=balance_triggered_global_unique_rematching_sensitivity",
    "primary_endpoint=n_sig_excluding_gKO_padj_0_05",
    "candidate_count=10",
    "controls_per_candidate=20",
    "control_reuse_across_candidates=FALSE",
    "matching_used_perturbation_endpoint=FALSE",
    paste0(
      "candidates_passing_individual_gate=",
      sum(candidates$individual_gate_pass)
    ),
    paste0("eligible_pooled_null_residuals=", nrow(finite_null)),
    paste0(
      "kruskal_wallis_p=",
      format(homogeneity$p.value, digits = 12)
    ),
    paste0("global_homogeneity_gate_pass=", global_gate),
    paste0(
      "candidates_with_empirical_p=",
      sum(is.finite(candidates$empirical_calibration_p))
    ),
    "TERT_status=not_calibration_eligible_due_to_absence_of_covariate_support",
    "claim_boundary=no_causal_or_target_specificity_claim",
    "original_v3_primary_calibration_unchanged=TRUE"
  ),
  gate_out,
  useBytes = TRUE
)

cat("GLOBAL_UNIQUE_CALIBRATION_SENSITIVITY_V1_COMPLETE\n")
cat(candidate_out, "\n", null_out, "\n", gate_out, "\n")
