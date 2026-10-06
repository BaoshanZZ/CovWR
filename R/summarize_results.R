summarize_simulation_results <- function(raw_results) {
  raw <- data.table::as.data.table(raw_results)

  scenario_cols <- intersect(
    c(
      "scenario_id", "N", "e", "effect", "prognosis",
      "outcome_shape", "effect_heterogeneity", "endpoint_prognosis_pattern",
      "n_noise_covariates", "adjustment_covariate_pattern", "scenario_change",
      "specification"
    ),
    names(raw)
  )
  metric_cols <- c(scenario_cols, "method", "estimand")

  has_bootstrap <- all(c(
    "bootstrap_B", "boot_var", "se_boot", "covered_boot", "reject_boot"
  ) %in% names(raw))

  summary <- raw[, .(
    truth = unique(truth)[1L],
    estimate_mean = mean(estimate, na.rm = TRUE),
    bias = mean(estimate - truth, na.rm = TRUE),
    ESD = stats::sd(estimate, na.rm = TRUE),
    ASE = mean(se, na.rm = TRUE),
    RMSE_SE = sqrt(mean(se^2, na.rm = TRUE)),
    CP = mean(covered, na.rm = TRUE),
    reject_rate = mean(reject, na.rm = TRUE),
    n_rep = .N
  ), by = metric_cols]

  if (has_bootstrap) {
    boot_summary <- raw[, .(
      bootstrap_B = max(bootstrap_B, na.rm = TRUE),
      Mean_Boot_Var = mean(boot_var, na.rm = TRUE),
      ASE_boot = mean(se_boot, na.rm = TRUE),
      RMSE_SE_boot = sqrt(mean(se_boot^2, na.rm = TRUE)),
      CP_boot = mean(covered_boot, na.rm = TRUE),
      reject_rate_boot = mean(reject_boot, na.rm = TRUE)
    ), by = metric_cols]
    boot_summary[!is.finite(bootstrap_B), bootstrap_B := 0L]
    boot_summary[!is.finite(Mean_Boot_Var), Mean_Boot_Var := NA_real_]
    boot_summary[!is.finite(ASE_boot), ASE_boot := NA_real_]
    boot_summary[!is.finite(RMSE_SE_boot), RMSE_SE_boot := NA_real_]
    boot_summary[!is.finite(CP_boot), CP_boot := NA_real_]
    boot_summary[!is.finite(reject_rate_boot), reject_rate_boot := NA_real_]
    summary <- merge(summary, boot_summary, by = metric_cols, all.x = TRUE, sort = FALSE)
  } else {
    summary[, `:=`(
      bootstrap_B = 0L,
      Mean_Boot_Var = NA_real_,
      ASE_boot = NA_real_,
      RMSE_SE_boot = NA_real_,
      CP_boot = NA_real_,
      reject_rate_boot = NA_real_
    )]
  }

  summary[, `:=`(
    MC_SE_mean = ESD / sqrt(n_rep),
    ASE_ESD = ASE / ESD
  )]

  unadj <- summary[method == "UNADJ", .(
    scenario_id,
    estimand,
    ESD_UNADJ = ESD
  )]
  summary <- merge(summary, unadj, by = c("scenario_id", "estimand"), all.x = TRUE)
  summary[, PVR := 1 - (ESD^2 / ESD_UNADJ^2)]
  summary[method == "UNADJ", PVR := 0]
  summary[, PVR_percent := 100 * PVR]
  summary[, ESD_UNADJ := NULL]

  data.table::setorder(summary, scenario_id, estimand, method)
  summary[]
}

check_simulation_outputs <- function(raw_results, summary_results) {
  forbidden <- c("IPW", "EB", "Entropy", "Flexible", "Proposed--Flexible", "Proposed Linear")
  bad_methods <- unique(raw_results$method[raw_results$method %in% forbidden])
  if (length(bad_methods) > 0L) {
    stop("Unexpected skipped method(s) found: ", paste(bad_methods, collapse = ", "), call. = FALSE)
  }
  allowed_methods <- c("UNADJ", "LADJ", "FADJ")
  bad_allowed <- setdiff(unique(raw_results$method), allowed_methods)
  if (length(bad_allowed) > 0L) {
    stop("Unexpected comparison method(s) found: ", paste(bad_allowed, collapse = ", "), call. = FALSE)
  }

  required_summary_cols <- c(
    "scenario_id", "N", "e", "effect", "prognosis",
    "method", "estimand", "truth", "estimate_mean", "bias", "ESD",
    "MC_SE_mean", "ASE", "ASE_ESD", "RMSE_SE", "CP", "reject_rate",
    "PVR", "PVR_percent"
  )
  if ("outcome_shape" %in% names(summary_results)) {
    required_summary_cols <- c(required_summary_cols, "outcome_shape")
  }
  if ("effect_heterogeneity" %in% names(summary_results)) {
    required_summary_cols <- c(required_summary_cols, "effect_heterogeneity")
  }
  if ("endpoint_prognosis_pattern" %in% names(summary_results)) {
    required_summary_cols <- c(required_summary_cols, "endpoint_prognosis_pattern")
  }
  if ("n_noise_covariates" %in% names(summary_results)) {
    required_summary_cols <- c(required_summary_cols, "n_noise_covariates")
  }
  if ("adjustment_covariate_pattern" %in% names(summary_results)) {
    required_summary_cols <- c(required_summary_cols, "adjustment_covariate_pattern")
  }
  missing_cols <- setdiff(required_summary_cols, names(summary_results))
  if (length(missing_cols) > 0L) {
    stop("Summary is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (anyNA(summary_results[, ..required_summary_cols])) {
    stop("Summary contains missing values in required columns.", call. = FALSE)
  }

  optional_boot_cols <- c(
    "bootstrap_B", "Mean_Boot_Var", "ASE_boot", "RMSE_SE_boot",
    "CP_boot", "reject_rate_boot"
  )
  missing_boot_cols <- setdiff(optional_boot_cols, names(summary_results))
  if (length(missing_boot_cols) > 0L) {
    stop("Summary is missing bootstrap column(s): ", paste(missing_boot_cols, collapse = ", "), call. = FALSE)
  }

  boot_rows <- summary_results$bootstrap_B > 0L
  if (any(boot_rows) && anyNA(summary_results[boot_rows, ..optional_boot_cols])) {
    stop("Summary contains missing values in bootstrap columns for bootstrap runs.", call. = FALSE)
  }

  invisible(TRUE)
}
