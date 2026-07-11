summarize_simulation_results <- function(raw_results) {
  raw <- data.table::as.data.table(raw_results)

  scenario_cols <- intersect(
    c(
      "scenario_id", "N", "e", "effect", "prognosis",
      "outcome_shape", "effect_heterogeneity", "endpoint_prognosis_pattern", "specification"
    ),
    names(raw)
  )
  metric_cols <- c(scenario_cols, "method", "estimand")

  summary <- raw[, .(
    truth = unique(truth)[1L],
    estimate_mean = mean(estimate, na.rm = TRUE),
    bias = mean(estimate - truth, na.rm = TRUE),
    ESD = stats::sd(estimate, na.rm = TRUE),
    ASE = mean(se, na.rm = TRUE),
    CP = mean(covered, na.rm = TRUE),
    reject_rate = mean(reject, na.rm = TRUE),
    n_rep = .N
  ), by = metric_cols]
  summary[, `:=`(
    MC_SE_mean = ESD / sqrt(n_rep)
  )]

  unadj <- summary[method == "UNADJ", .(
    scenario_id,
    estimand,
    ESD_UNADJ = ESD
  )]
  summary <- merge(summary, unadj, by = c("scenario_id", "estimand"), all.x = TRUE)
  summary[, PVR := 1 - (ESD^2 / ESD_UNADJ^2)]
  summary[method == "UNADJ", PVR := 0]
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
    "MC_SE_mean", "ASE", "CP", "reject_rate", "PVR"
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
  missing_cols <- setdiff(required_summary_cols, names(summary_results))
  if (length(missing_cols) > 0L) {
    stop("Summary is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (anyNA(summary_results[, ..required_summary_cols])) {
    stop("Summary contains missing values in required columns.", call. = FALSE)
  }

  invisible(TRUE)
}
