options(stringsAsFactors = FALSE, scipen = 999)
Sys.setenv(TZ = "America/New_York")

known_e <- 0.5
analysis_horizon <- 28
primary_covariates <- c(
  "X_SOFA_0d",
  "X_SIRS_0d",
  "X_Ddimer_0d",
  "X_platelet_0d"
)
sensitivity_label <- paste(
  "LADJ sensitivity including baseline TM with pooled imputation",
  "+ missing indicator"
)

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", file_arg)))
  }
  normalizePath("analysis_trial23_win_ladj.R")
}

analysis_dir <- dirname(script_path())
data_path <- file.path(analysis_dir, "trial23.rds")
dictionary_path <- file.path(analysis_dir, "trial23_dictionary.csv")
output_dir <- file.path(analysis_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_analysis <- function(...) stop(..., call. = FALSE)
assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop_analysis(message)
}

format_number <- function(x, digits = 4L) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
}

markdown_table <- function(x, digits = 4L) {
  x <- as.data.frame(x, check.names = FALSE)
  if (nrow(x) == 0L) return("_No rows._")
  formatted <- lapply(x, function(column) {
    if (is.numeric(column)) format_number(column, digits) else as.character(column)
  })
  formatted <- as.data.frame(formatted, check.names = FALSE)
  formatted[is.na(formatted)] <- ""
  header <- paste0("| ", paste(names(formatted), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(formatted)), collapse = " | "), " |")
  rows <- apply(formatted, 1L, function(row) {
    paste0("| ", paste(gsub("\\|", "\\\\|", row), collapse = " | "), " |")
  })
  paste(c(header, separator, rows), collapse = "\n")
}

write_csv <- function(x, filename) {
  utils::write.csv(x, file.path(output_dir, filename), row.names = FALSE, na = "")
}

assert_true(file.exists(data_path), paste("Missing RDS file:", data_path))
assert_true(file.exists(dictionary_path), paste("Missing dictionary:", dictionary_path))
dat <- as.data.frame(readRDS(data_path), check.names = FALSE)
dictionary <- utils::read.csv(dictionary_path, check.names = FALSE)
N <- nrow(dat)
dat$participant_id <- seq_len(N)

required_outcomes <- c(
  "YP_event_28d",
  "YP_time_28d",
  "YS_delta_SOFA_3d",
  "YS_delta_Ddimer_3d"
)
required_variables <- c("Treatment", required_outcomes, primary_covariates, "X_TM")
missing_variables <- setdiff(required_variables, names(dat))
assert_true(
  length(missing_variables) == 0L,
  paste("Required variable(s) missing:", paste(missing_variables, collapse = ", "))
)

expected_treatment_levels <- c("No rhTM", "rhTM")
treatment_levels <- unique(as.character(dat$Treatment))
assert_true(
  setequal(treatment_levels, expected_treatment_levels),
  paste("Unexpected treatment levels:", paste(treatment_levels, collapse = ", "))
)
dat$A <- as.integer(as.character(dat$Treatment) == "rhTM")
n1 <- sum(dat$A == 1L)
n0 <- sum(dat$A == 0L)
assert_true(n1 > 0L && n0 > 0L, "Both randomized arms must be present.")

assert_true(!anyNA(dat$YP_event_28d), "YP_event_28d contains missing values.")
assert_true(!anyNA(dat$YP_time_28d), "YP_time_28d contains missing values.")
event_character <- as.character(dat$YP_event_28d)
event_levels <- unique(event_character)
early_time_levels <- unique(event_character[dat$YP_time_28d < analysis_horizon])
assert_true(
  length(early_time_levels) == 1L,
  "Could not uniquely determine the death/event value from times before Day 28."
)
death_value <- early_time_levels[[1L]]
dat$death_event <- as.integer(event_character == death_value)
assert_true(
  all(dat$YP_time_28d[dat$death_event == 0L] == analysis_horizon),
  "A presumed Day-28 survivor has follow-up time shorter than 28 days."
)
assert_true(
  all(dat$YP_time_28d[dat$death_event == 1L] < analysis_horizon),
  "A presumed death/event does not occur before Day 28."
)
assert_true(
  all(!is.na(dat$YS_delta_SOFA_3d[dat$death_event == 0L])),
  "At least one Day-28 survivor is missing Day-3 SOFA; primary analysis stopped."
)
assert_true(
  all(vapply(dat[primary_covariates], function(x) !anyNA(x), logical(1))),
  "A primary LADJ baseline covariate contains missing values."
)

prefixed_variables <- grep("^(YP_|YS_|X_)", names(dat), value = TRUE)
audit_variables <- c("Treatment", prefixed_variables)
audit_row <- function(variable) {
  x <- dat[[variable]]
  nonmissing <- x[!is.na(x)]
  numeric_x <- is.numeric(x) || is.integer(x)
  data.frame(
    Variable = variable,
    Class = paste(class(x), collapse = ";"),
    N = length(x),
    Missing = sum(is.na(x)),
    Unique_nonmissing = length(unique(nonmissing)),
    Minimum = if (numeric_x && length(nonmissing)) min(nonmissing) else NA_real_,
    Maximum = if (numeric_x && length(nonmissing)) max(nonmissing) else NA_real_,
    Levels_or_values = if (numeric_x) {
      if (length(unique(nonmissing)) <= 20L) paste(sort(unique(nonmissing)), collapse = "; ") else ""
    } else {
      paste(unique(as.character(nonmissing)), collapse = "; ")
    },
    Used_in_primary = variable %in% c("Treatment", required_outcomes, primary_covariates),
    stringsAsFactors = FALSE
  )
}
data_audit <- do.call(rbind, lapply(audit_variables, audit_row))
dictionary_columns <- intersect(
  c("variable_name", "variable_role", "variable_type", "short_explanation"),
  names(dictionary)
)
if (length(dictionary_columns) > 1L) {
  dictionary_subset <- dictionary[dictionary_columns]
  data_audit <- merge(
    data_audit,
    dictionary_subset,
    by.x = "Variable",
    by.y = "variable_name",
    all.x = TRUE,
    sort = FALSE
  )
  data_audit <- data_audit[match(audit_variables, data_audit$Variable), , drop = FALSE]
}
write_csv(data_audit, "trial23_data_audit.csv")

cat("Trial 23 input audit\n")
cat("N =", N, "; n1 (rhTM) =", n1, "; n0 (No rhTM) =", n0, "\n")
cat("Treatment levels:", paste(levels(dat$Treatment), collapse = ", "), "\n")
cat("Known randomization probability e =", known_e, "\n")
cat("Death/event value of YP_event_28d =", death_value, "\n")
cat("YP/YS/X variables:\n  ", paste(prefixed_variables, collapse = "\n  "), "\n", sep = "")
cat("Missing counts for analysis variables:\n")
print(data_audit[data_audit$Variable %in% required_variables, c("Variable", "Missing")], row.names = FALSE)
cat("Outcome ranges and unique values:\n")
print(data_audit[data_audit$Variable %in% required_outcomes,
                 c("Variable", "Minimum", "Maximum", "Levels_or_values")], row.names = FALSE)
cat("Clinical direction: lower Day-3 delta SOFA and lower Day-3 delta D-dimer are favorable.\n\n")

make_pairwise_comparisons <- function(data, sofa_rule = c("threshold_free", "margin_2")) {
  sofa_rule <- match.arg(sofa_rule)
  treated_rows <- which(data$A == 1L)
  control_rows <- which(data$A == 0L)
  n_t <- length(treated_rows)
  n_c <- length(control_rows)
  dims <- c(n_t, n_c)
  win <- matrix(0L, nrow = n_t, ncol = n_c)
  loss <- matrix(0L, nrow = n_t, ncol = n_c)
  tie <- matrix(0L, nrow = n_t, ncol = n_c)
  decision_tier <- matrix(NA_character_, nrow = n_t, ncol = n_c)

  tr_death <- data$death_event[treated_rows] == 1L
  co_death <- data$death_event[control_rows] == 1L
  tr_time <- data$YP_time_28d[treated_rows]
  co_time <- data$YP_time_28d[control_rows]
  tr_death_mat <- matrix(tr_death, nrow = n_t, ncol = n_c)
  co_death_mat <- matrix(co_death, nrow = n_t, ncol = n_c, byrow = TRUE)
  tr_time_mat <- matrix(tr_time, nrow = n_t, ncol = n_c)
  co_time_mat <- matrix(co_time, nrow = n_t, ncol = n_c, byrow = TRUE)

  survival_win <- !tr_death_mat & co_death_mat
  survival_loss <- tr_death_mat & !co_death_mat
  both_die <- tr_death_mat & co_death_mat
  survival_win <- survival_win | (both_die & tr_time_mat > co_time_mat)
  survival_loss <- survival_loss | (both_die & tr_time_mat < co_time_mat)
  same_death_time <- both_die & tr_time_mat == co_time_mat
  both_survive <- !tr_death_mat & !co_death_mat

  win[survival_win] <- 1L
  loss[survival_loss] <- 1L
  tie[same_death_time] <- 1L
  decision_tier[survival_win | survival_loss] <- "28-day survival"
  decision_tier[same_death_time] <- "Final tie at survival"

  tr_sofa <- data$YS_delta_SOFA_3d[treated_rows]
  co_sofa <- data$YS_delta_SOFA_3d[control_rows]
  tr_sofa_mat <- matrix(tr_sofa, nrow = n_t, ncol = n_c)
  co_sofa_mat <- matrix(co_sofa, nrow = n_t, ncol = n_c, byrow = TRUE)
  assert_true(
    !anyNA(tr_sofa_mat[both_survive]) && !anyNA(co_sofa_mat[both_survive]),
    paste("Missing Day-3 SOFA in a surviving pair under", sofa_rule)
  )
  if (sofa_rule == "threshold_free") {
    sofa_win <- both_survive & tr_sofa_mat < co_sofa_mat
    sofa_loss <- both_survive & tr_sofa_mat > co_sofa_mat
    pass_sofa <- both_survive & tr_sofa_mat == co_sofa_mat
    definition_label <- "Threshold-free SOFA"
  } else {
    sofa_win <- both_survive & tr_sofa_mat <= co_sofa_mat - 2
    sofa_loss <- both_survive & co_sofa_mat <= tr_sofa_mat - 2
    pass_sofa <- both_survive & abs(tr_sofa_mat - co_sofa_mat) < 2
    definition_label <- "2-point SOFA margin"
  }
  win[sofa_win] <- 1L
  loss[sofa_loss] <- 1L
  decision_tier[sofa_win | sofa_loss] <- "Day-3 SOFA"

  tr_ddimer <- data$YS_delta_Ddimer_3d[treated_rows]
  co_ddimer <- data$YS_delta_Ddimer_3d[control_rows]
  tr_ddimer_mat <- matrix(tr_ddimer, nrow = n_t, ncol = n_c)
  co_ddimer_mat <- matrix(co_ddimer, nrow = n_t, ncol = n_c, byrow = TRUE)
  missing_ddimer_pair <- pass_sofa & (is.na(tr_ddimer_mat) | is.na(co_ddimer_mat))
  assert_true(
    !any(missing_ddimer_pair),
    paste("A pair reaching Day-3 D-dimer has missing D-dimer under", definition_label)
  )
  ddimer_win <- pass_sofa & tr_ddimer_mat < co_ddimer_mat
  ddimer_loss <- pass_sofa & tr_ddimer_mat > co_ddimer_mat
  ddimer_tie <- pass_sofa & tr_ddimer_mat == co_ddimer_mat
  win[ddimer_win] <- 1L
  loss[ddimer_loss] <- 1L
  tie[ddimer_tie] <- 1L
  decision_tier[ddimer_win | ddimer_loss] <- "Day-3 D-dimer"
  decision_tier[ddimer_tie] <- "Final tie at D-dimer"

  assert_true(all(win + loss + tie == 1L), paste("Pair partition failed for", definition_label))
  assert_true(!anyNA(decision_tier), paste("Decision tier missing for", definition_label))

  pair_index <- expand.grid(
    treated_row = treated_rows,
    control_row = control_rows,
    KEEP.OUT.ATTRS = FALSE
  )
  pairs <- data.frame(
    HCE_definition = definition_label,
    treated_id = data$participant_id[pair_index$treated_row],
    control_id = data$participant_id[pair_index$control_row],
    win_ij = as.integer(win),
    loss_ij = as.integer(loss),
    tie_ij = as.integer(tie),
    decision_tier = as.character(decision_tier),
    stringsAsFactors = FALSE
  )
  list(
    label = definition_label,
    win = win,
    loss = loss,
    tie = tie,
    decision_tier = decision_tier,
    pairs = pairs,
    treated_rows = treated_rows,
    control_rows = control_rows,
    pairs_reaching_ddimer = sum(pass_sofa),
    sofa_decided = sofa_win | sofa_loss
  )
}

measure_gradient <- function(measure, tau_w, tau_l) {
  delta <- tau_w - tau_l
  switch(
    measure,
    NB = c(1, -1),
    WR = c(1 / tau_w, -1 / tau_l),
    WO = 2 / (1 - delta^2) * c(1, -1),
    stop_analysis("Unknown measure: ", measure)
  )
}

measure_values <- function(tau_w, tau_l) {
  assert_true(tau_w > 0 && tau_l > 0, "Win/loss probability must be positive for WR.")
  delta <- tau_w - tau_l
  assert_true(abs(delta) < 1, "Absolute net benefit must be less than one for WO.")
  c(
    NB = delta,
    WR = log(tau_w) - log(tau_l),
    WO = log1p(delta) - log1p(-delta)
  )
}

make_inference_rows <- function(hce_label, method, tau, psi_components) {
  values <- measure_values(tau[["w"]], tau[["l"]])
  psi_components <- as.matrix(psi_components)
  centered <- sweep(psi_components, 2L, colMeans(psi_components), "-")
  sigma <- crossprod(centered) / N
  measures <- names(values)
  rows <- lapply(measures, function(measure) {
    gradient <- measure_gradient(measure, tau[["w"]], tau[["l"]])
    variance <- as.numeric(crossprod(gradient, sigma %*% gradient) / N)
    assert_true(is.finite(variance) && variance >= -1e-14, "Invalid estimated variance.")
    variance <- max(variance, 0)
    se <- sqrt(variance)
    estimate_log <- unname(values[[measure]])
    ci_log <- estimate_log + stats::qnorm(c(0.025, 0.975)) * se
    p_value <- 2 * stats::pnorm(-abs(estimate_log / se))
    if (measure == "NB") {
      estimate_report <- estimate_log
      ci_report <- ci_log
      estimate_log_output <- NA_real_
    } else {
      estimate_report <- exp(estimate_log)
      ci_report <- exp(ci_log)
      estimate_log_output <- estimate_log
    }
    data.frame(
      HCE_definition = hce_label,
      Method = method,
      Measure = measure,
      Estimate_log_scale = estimate_log_output,
      Estimate_report_scale = estimate_report,
      SE = se,
      CI_lower = ci_report[[1L]],
      CI_upper = ci_report[[2L]],
      CI_lower_log_scale = if (measure == "NB") NA_real_ else ci_log[[1L]],
      CI_upper_log_scale = if (measure == "NB") NA_real_ else ci_log[[2L]],
      P_value = p_value,
      Estimated_variance = variance,
      tau_w = tau[["w"]],
      tau_l = tau[["l"]],
      tau_tie = 1 - tau[["w"]] - tau[["l"]],
      stringsAsFactors = FALSE
    )
  })
  list(rows = do.call(rbind, rows), sigma = sigma)
}

construct_scores <- function(pair_object) {
  W <- pair_object$win
  L <- pair_object$loss
  tau_w <- mean(W)
  tau_l <- mean(L)
  placement <- list(
    g_w_T = rowMeans(W) - tau_w,
    g_l_T = rowMeans(L) - tau_l,
    g_w_C = colMeans(W) - tau_w,
    g_l_C = colMeans(L) - tau_l
  )
  psi_w <- numeric(N)
  psi_l <- numeric(N)
  psi_w[pair_object$treated_rows] <- placement$g_w_T / known_e
  psi_l[pair_object$treated_rows] <- placement$g_l_T / known_e
  psi_w[pair_object$control_rows] <- placement$g_w_C / (1 - known_e)
  psi_l[pair_object$control_rows] <- placement$g_l_C / (1 - known_e)
  list(
    tau = c(w = tau_w, l = tau_l),
    placement = placement,
    psi = cbind(w = psi_w, l = psi_l)
  )
}

fit_ladj <- function(scores, X, covariate_names, method_label = "LADJ") {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  assert_true(!anyNA(X), paste("Missing LADJ feature in", method_label))
  X_center <- sweep(X, 2L, colMeans(X), "-")
  colnames(X_center) <- covariate_names
  Z <- dat$A - known_e
  design <- X_center * Z
  gram <- crossprod(design)
  rank <- qr(gram)$rank
  eigenvalues <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
  smallest_eigenvalue <- min(eigenvalues)
  condition_number <- max(eigenvalues) / smallest_eigenvalue
  assert_true(rank == ncol(X_center), paste("Singular LADJ Gram matrix in", method_label))
  assert_true(
    smallest_eigenvalue > .Machine$double.eps * max(eigenvalues) * ncol(X_center),
    paste("Numerically singular LADJ Gram matrix in", method_label)
  )
  beta <- solve(gram, crossprod(design, scores$psi))
  colnames(beta) <- c("beta_w", "beta_l")
  fitted_b <- X_center %*% beta
  psi_adjusted <- scores$psi - fitted_b * Z
  direct_correction <- colMeans(fitted_b * Z)

  mean_treated <- colMeans(X[dat$A == 1L, , drop = FALSE])
  mean_control <- colMeans(X[dat$A == 0L, , drop = FALSE])
  imbalance <- mean_treated - mean_control
  imbalance_factor <- n1 * n0 / N^2
  identity_left <- colMeans(X_center * Z)
  identity_right <- imbalance_factor * imbalance
  imbalance_correction <- as.vector(imbalance_factor * crossprod(imbalance, beta))
  tau_adjusted <- scores$tau - direct_correction

  feature_sd <- apply(X, 2L, stats::sd)
  assert_true(all(feature_sd > 0), paste("Zero-variance LADJ feature in", method_label))
  X_scaled <- sweep(X_center, 2L, feature_sd, "/")
  design_scaled <- X_scaled * Z
  gram_scaled <- crossprod(design_scaled)
  beta_scaled <- solve(gram_scaled, crossprod(design_scaled, scores$psi))
  fitted_scaled <- X_scaled %*% beta_scaled
  scaled_eigenvalues <- eigen(gram_scaled, symmetric = TRUE, only.values = TRUE)$values

  diagnostics <- data.frame(
    Method = method_label,
    Number_features = ncol(X_center),
    Gram_rank = rank,
    Smallest_eigenvalue = smallest_eigenvalue,
    Largest_eigenvalue = max(eigenvalues),
    Condition_number = condition_number,
    Scaled_condition_number = max(scaled_eigenvalues) / min(scaled_eigenvalues),
    Max_centering_identity_error = max(abs(identity_left - identity_right)),
    Max_tau_formula_error = max(abs(direct_correction - imbalance_correction)),
    Max_scaled_fitted_value_error = max(abs(fitted_b - fitted_scaled)),
    stringsAsFactors = FALSE
  )
  coefficients <- data.frame(
    Method = method_label,
    Covariate = covariate_names,
    beta_w = beta[, "beta_w"],
    beta_l = beta[, "beta_l"],
    stringsAsFactors = FALSE
  )
  list(
    tau = tau_adjusted,
    psi = psi_adjusted,
    beta = beta,
    fitted_b = fitted_b,
    coefficients = coefficients,
    diagnostics = diagnostics,
    covariance = crossprod(sweep(psi_adjusted, 2L, colMeans(psi_adjusted), "-")) / N
  )
}

primary_X <- dat[primary_covariates]
tm_median <- stats::median(dat$X_TM, na.rm = TRUE)
sensitivity_X <- cbind(
  primary_X,
  X_TM_imp = ifelse(is.na(dat$X_TM), tm_median, dat$X_TM),
  X_TM_missing = as.integer(is.na(dat$X_TM))
)

pooled_within_arm_sd <- function(x, A) {
  n_t <- sum(A == 1L)
  n_c <- sum(A == 0L)
  sqrt(((n_t - 1) * stats::var(x[A == 1L]) + (n_c - 1) * stats::var(x[A == 0L])) /
         (n_t + n_c - 2))
}
baseline_balance <- do.call(rbind, lapply(primary_covariates, function(variable) {
  x <- dat[[variable]]
  mean_t <- mean(x[dat$A == 1L])
  mean_c <- mean(x[dat$A == 0L])
  pooled_sd <- pooled_within_arm_sd(x, dat$A)
  data.frame(
    Covariate = variable,
    Treatment_mean = mean_t,
    Control_mean = mean_c,
    Raw_mean_difference = mean_t - mean_c,
    Pooled_within_arm_SD = pooled_sd,
    Overall_SD = stats::sd(x),
    Standardized_mean_difference = (mean_t - mean_c) / pooled_sd,
    stringsAsFactors = FALSE
  )
}))
write_csv(baseline_balance, "trial23_baseline_balance.csv")

pair_objects <- list(
  threshold_free = make_pairwise_comparisons(dat, "threshold_free"),
  margin_2 = make_pairwise_comparisons(dat, "margin_2")
)
reassigned_to_ddimer <- sum(
  pair_objects$threshold_free$sofa_decided &
    !pair_objects$margin_2$sofa_decided
)

all_pair_rows <- list()
pairwise_summary_rows <- list()
all_result_rows <- list()
all_efficiency_rows <- list()
all_coefficient_rows <- list()
all_diagnostic_rows <- list()
all_score_rows <- list()
all_qc_rows <- list()

add_qc <- function(hce, check, observed, tolerance, passed, details = "") {
  data.frame(
    HCE_definition = hce,
    Check = check,
    Observed = observed,
    Tolerance = tolerance,
    Passed = passed,
    Details = details,
    stringsAsFactors = FALSE
  )
}

for (object_name in names(pair_objects)) {
  pair_object <- pair_objects[[object_name]]
  hce <- pair_object$label
  scores <- construct_scores(pair_object)
  primary_fit <- fit_ladj(scores, primary_X, primary_covariates, "LADJ")
  sensitivity_fit <- fit_ladj(
    scores,
    sensitivity_X,
    colnames(sensitivity_X),
    sensitivity_label
  )

  unadjusted_inference <- make_inference_rows(hce, "UNADJ", scores$tau, scores$psi)
  ladj_inference <- make_inference_rows(hce, "LADJ", primary_fit$tau, primary_fit$psi)
  sensitivity_inference <- make_inference_rows(
    hce,
    sensitivity_label,
    sensitivity_fit$tau,
    sensitivity_fit$psi
  )
  method_results <- rbind(
    unadjusted_inference$rows,
    ladj_inference$rows,
    sensitivity_inference$rows
  )

  unadjusted_variance <- setNames(
    unadjusted_inference$rows$Estimated_variance,
    unadjusted_inference$rows$Measure
  )
  efficiency_rows <- lapply(c("LADJ", sensitivity_label), function(method) {
    method_rows <- method_results[method_results$Method == method, , drop = FALSE]
    method_variance <- setNames(method_rows$Estimated_variance, method_rows$Measure)
    data.frame(
      HCE_definition = hce,
      Method = method,
      Measure = names(unadjusted_variance),
      Unadjusted_variance = as.numeric(unadjusted_variance),
      Adjusted_variance = as.numeric(method_variance[names(unadjusted_variance)]),
      Variance_reduction_percent = 100 * (1 - method_variance[names(unadjusted_variance)] /
                                            unadjusted_variance),
      Relative_efficiency = unadjusted_variance /
        method_variance[names(unadjusted_variance)],
      stringsAsFactors = FALSE
    )
  })
  efficiency <- do.call(rbind, efficiency_rows)
  method_results$Variance_reduction_percent <- 0
  method_results$Relative_efficiency <- 1
  for (row_id in seq_len(nrow(efficiency))) {
    idx <- method_results$Method == efficiency$Method[[row_id]] &
      method_results$Measure == efficiency$Measure[[row_id]]
    method_results$Variance_reduction_percent[idx] <- efficiency$Variance_reduction_percent[[row_id]]
    method_results$Relative_efficiency[idx] <- efficiency$Relative_efficiency[[row_id]]
  }

  pairs <- pair_object$pairs
  N_pairs <- nrow(pairs)
  decided_survival <- sum(pairs$decision_tier == "28-day survival")
  decided_sofa <- sum(pairs$decision_tier == "Day-3 SOFA")
  decided_ddimer <- sum(pairs$decision_tier == "Day-3 D-dimer")
  final_ties <- sum(pairs$tie_ij)
  ties_survival <- sum(pairs$decision_tier == "Final tie at survival")
  ties_ddimer <- sum(pairs$decision_tier == "Final tie at D-dimer")
  reassigned <- if (object_name == "margin_2") reassigned_to_ddimer else 0L
  pair_summary <- data.frame(
    HCE_definition = hce,
    N = N,
    n1 = n1,
    n0 = n0,
    N_pairs = N_pairs,
    Treatment_wins = sum(pairs$win_ij),
    Treatment_losses = sum(pairs$loss_ij),
    Final_ties = final_ties,
    tau_w = mean(pairs$win_ij),
    tau_l = mean(pairs$loss_ij),
    tau_tie = mean(pairs$tie_ij),
    Decided_by_survival = decided_survival,
    Decided_by_survival_percent = 100 * decided_survival / N_pairs,
    Decided_by_SOFA = decided_sofa,
    Decided_by_SOFA_percent = 100 * decided_sofa / N_pairs,
    Decided_by_Ddimer = decided_ddimer,
    Decided_by_Ddimer_percent = 100 * decided_ddimer / N_pairs,
    Final_ties_percent = 100 * final_ties / N_pairs,
    Final_ties_at_survival = ties_survival,
    Final_ties_at_Ddimer = ties_ddimer,
    Reached_Ddimer_tier = pair_object$pairs_reaching_ddimer,
    SOFA_pairs_reassigned_to_Ddimer_vs_threshold_free = reassigned,
    SOFA_pairs_reassigned_percent_of_all_pairs = 100 * reassigned / N_pairs,
    stringsAsFactors = FALSE
  )

  primary_fit$coefficients$HCE_definition <- hce
  sensitivity_fit$coefficients$HCE_definition <- hce
  primary_fit$diagnostics$HCE_definition <- hce
  sensitivity_fit$diagnostics$HCE_definition <- hce
  score_rows <- data.frame(
    HCE_definition = hce,
    participant_id = dat$participant_id,
    Treatment = as.character(dat$Treatment),
    A = dat$A,
    psi_w = scores$psi[, "w"],
    psi_l = scores$psi[, "l"],
    psi_w_LADJ = primary_fit$psi[, "w"],
    psi_l_LADJ = primary_fit$psi[, "l"],
    stringsAsFactors = FALSE
  )

  covariance_eigen <- eigen(primary_fit$covariance, symmetric = TRUE, only.values = TRUE)$values
  placement_means <- c(
    mean(scores$placement$g_w_T),
    mean(scores$placement$g_l_T),
    mean(scores$placement$g_w_C),
    mean(scores$placement$g_l_C)
  )
  qc <- rbind(
    add_qc(hce, "Number of pairs equals n1*n0", N_pairs, n1 * n0, N_pairs == n1 * n0),
    add_qc(hce, "Each pair has win+loss+tie=1", max(abs(pairs$win_ij + pairs$loss_ij + pairs$tie_ij - 1)), 0, all(pairs$win_ij + pairs$loss_ij + pairs$tie_ij == 1)),
    add_qc(hce, "Pair probabilities sum to one", abs(sum(scores$tau) + mean(pairs$tie_ij) - 1), 1e-12, abs(sum(scores$tau) + mean(pairs$tie_ij) - 1) < 1e-12),
    add_qc(hce, "Placement means are zero", max(abs(placement_means)), 1e-12, max(abs(placement_means)) < 1e-12),
    add_qc(hce, "Mean elementary IFs are zero", max(abs(colMeans(scores$psi))), 1e-12, max(abs(colMeans(scores$psi))) < 1e-12),
    add_qc(hce, "Known e is exactly 0.5", known_e, 0.5, identical(known_e, 0.5)),
    add_qc(hce, "No intercept in primary LADJ", as.numeric("(Intercept)" %in% primary_covariates), 0, !("(Intercept)" %in% primary_covariates)),
    add_qc(hce, "Primary LADJ uses exactly four prespecified covariates", length(primary_covariates), 4, identical(colnames(primary_X), primary_covariates)),
    add_qc(hce, "Equivalent adjusted tau formulas agree", primary_fit$diagnostics$Max_tau_formula_error, 1e-12, primary_fit$diagnostics$Max_tau_formula_error < 1e-12),
    add_qc(hce, "Primary Gram matrix is full rank", primary_fit$diagnostics$Gram_rank, length(primary_covariates), primary_fit$diagnostics$Gram_rank == length(primary_covariates)),
    add_qc(hce, "Adjusted covariance is symmetric", max(abs(primary_fit$covariance - t(primary_fit$covariance))), 1e-12, max(abs(primary_fit$covariance - t(primary_fit$covariance))) < 1e-12),
    add_qc(hce, "Adjusted covariance is positive semidefinite", min(covariance_eigen), -1e-12, min(covariance_eigen) >= -1e-12),
    add_qc(hce, "No survivor missing Day-3 SOFA", sum(is.na(dat$YS_delta_SOFA_3d[dat$death_event == 0L])), 0, !anyNA(dat$YS_delta_SOFA_3d[dat$death_event == 0L])),
    add_qc(hce, "No pair reaching D-dimer has missing D-dimer", 0, 0, TRUE),
    add_qc(hce, "All randomized subjects retained", length(unique(c(pairs$treated_id, pairs$control_id))), N, length(unique(c(pairs$treated_id, pairs$control_id))) == N)
  )

  all_pair_rows[[hce]] <- pairs
  pairwise_summary_rows[[hce]] <- pair_summary
  all_result_rows[[hce]] <- method_results
  all_efficiency_rows[[hce]] <- efficiency
  all_coefficient_rows[[hce]] <- rbind(primary_fit$coefficients, sensitivity_fit$coefficients)
  all_diagnostic_rows[[hce]] <- rbind(primary_fit$diagnostics, sensitivity_fit$diagnostics)
  all_score_rows[[hce]] <- score_rows
  all_qc_rows[[hce]] <- qc
}

pair_comparisons <- do.call(rbind, all_pair_rows)
row.names(pair_comparisons) <- NULL
pairwise_summary <- do.call(rbind, pairwise_summary_rows)
row.names(pairwise_summary) <- NULL
results <- do.call(rbind, all_result_rows)
row.names(results) <- NULL
efficiency <- do.call(rbind, all_efficiency_rows)
row.names(efficiency) <- NULL
coefficients <- do.call(rbind, all_coefficient_rows)
row.names(coefficients) <- NULL
coefficients <- coefficients[c("HCE_definition", "Method", "Covariate", "beta_w", "beta_l")]
diagnostics <- do.call(rbind, all_diagnostic_rows)
row.names(diagnostics) <- NULL
diagnostics <- diagnostics[c("HCE_definition", setdiff(names(diagnostics), "HCE_definition"))]
subject_scores <- do.call(rbind, all_score_rows)
row.names(subject_scores) <- NULL
quality_control <- do.call(rbind, all_qc_rows)
row.names(quality_control) <- NULL

same_subjects <- identical(
  sort(unique(c(
    pair_objects$threshold_free$pairs$treated_id,
    pair_objects$threshold_free$pairs$control_id
  ))),
  sort(unique(c(
    pair_objects$margin_2$pairs$treated_id,
    pair_objects$margin_2$pairs$control_id
  )))
)
quality_control <- rbind(
  quality_control,
  add_qc(
    "Across HCE definitions",
    "Same subjects and baseline procedure; only SOFA rule changes",
    as.numeric(same_subjects),
    1,
    same_subjects && identical(primary_covariates, colnames(primary_X))
  )
)
assert_true(all(quality_control$Passed), "At least one quality-control check failed.")

primary_results <- results[results$Method %in% c("UNADJ", "LADJ"), ]
tm_sensitivity_results <- results[results$Method == sensitivity_label, ]
primary_efficiency <- efficiency[efficiency$Method == "LADJ", ]
tm_sensitivity_efficiency <- efficiency[efficiency$Method == sensitivity_label, ]
primary_coefficients <- coefficients[coefficients$Method == "LADJ", ]
tm_sensitivity_coefficients <- coefficients[coefficients$Method == sensitivity_label, ]
primary_diagnostics <- diagnostics[diagnostics$Method == "LADJ", ]
tm_sensitivity_diagnostics <- diagnostics[diagnostics$Method == sensitivity_label, ]

write_csv(pair_comparisons, "trial23_pairwise_comparisons.csv")
write_csv(pairwise_summary, "trial23_pairwise_summary.csv")
write_csv(primary_results, "trial23_unadjusted_ladj_results.csv")
write_csv(primary_efficiency, "trial23_efficiency_comparison.csv")
write_csv(primary_coefficients, "trial23_ladj_coefficients.csv")
write_csv(primary_diagnostics, "trial23_ladj_diagnostics.csv")
write_csv(tm_sensitivity_results, "trial23_tm_sensitivity_results.csv")
write_csv(tm_sensitivity_efficiency, "trial23_tm_sensitivity_efficiency.csv")
write_csv(tm_sensitivity_coefficients, "trial23_tm_sensitivity_coefficients.csv")
write_csv(tm_sensitivity_diagnostics, "trial23_tm_sensitivity_diagnostics.csv")
write_csv(subject_scores, "trial23_subject_scores.csv")
write_csv(quality_control, "trial23_quality_control.csv")

open_plot_device <- function(filename, width, height, res) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename, width = width, height = height, units = "px", res = res)
    return(invisible(filename))
  }
  pdf_filename <- sub("\\.png$", ".pdf", filename)
  grDevices::pdf(pdf_filename, width = width / res, height = height / res)
  warning("ragg is unavailable; wrote PDF instead of PNG: ", pdf_filename)
  invisible(pdf_filename)
}

open_plot_device(
  file.path(output_dir, "trial23_primary_effect_estimates.png"),
  width = 2100,
  height = 850,
  res = 180
)
graphics::par(mfrow = c(1, 3), mar = c(5, 12, 4, 2), oma = c(0, 0, 2, 0))
forest_labels <- paste(
  ifelse(primary_results$HCE_definition == "Threshold-free SOFA", "Threshold-free", "2-point margin"),
  primary_results$Method,
  sep = ": "
)
method_colors <- c(UNADJ = "#5B6770", LADJ = "#1F77B4")
for (measure in c("NB", "WR", "WO")) {
  plot_data <- primary_results[primary_results$Measure == measure, ]
  y <- rev(seq_len(nrow(plot_data)))
  xlim <- range(c(plot_data$CI_lower, plot_data$CI_upper))
  padding <- diff(xlim) * 0.08
  xlim <- xlim + c(-padding, padding)
  graphics::plot(
    plot_data$Estimate_report_scale,
    y,
    xlim = xlim,
    ylim = c(0.5, nrow(plot_data) + 0.5),
    log = if (measure == "NB") "" else "x",
    xlab = if (measure == "NB") "Net benefit" else measure,
    ylab = "",
    yaxt = "n",
    pch = 19,
    col = method_colors[plot_data$Method],
    main = measure
  )
  graphics::abline(v = if (measure == "NB") 0 else 1, lty = 2, col = "grey45")
  graphics::segments(
    plot_data$CI_lower,
    y,
    plot_data$CI_upper,
    y,
    lwd = 2,
    col = method_colors[plot_data$Method]
  )
  graphics::points(
    plot_data$Estimate_report_scale,
    y,
    pch = 19,
    cex = 1.1,
    col = method_colors[plot_data$Method]
  )
  graphics::axis(2, at = y, labels = forest_labels[primary_results$Measure == measure], las = 1)
}
graphics::mtext("Trial 23 primary estimates with 95% Wald confidence intervals", outer = TRUE, cex = 1.2)
grDevices::dev.off()

tier_matrix <- rbind(
  "28-day survival" = pairwise_summary$Decided_by_survival_percent,
  "Day-3 SOFA" = pairwise_summary$Decided_by_SOFA_percent,
  "Day-3 D-dimer" = pairwise_summary$Decided_by_Ddimer_percent,
  "Final tie" = pairwise_summary$Final_ties_percent
)
open_plot_device(
  file.path(output_dir, "trial23_pair_decision_tiers.png"),
  width = 1400,
  height = 900,
  res = 180
)
graphics::par(mar = c(6, 5, 4, 11), xpd = TRUE)
graphics::barplot(
  tier_matrix,
  names.arg = c("Threshold-free SOFA", "2-point SOFA margin"),
  col = c("#4C78A8", "#72B7B2", "#F2CF5B", "#B0B0B0"),
  ylim = c(0, 100),
  ylab = "Percentage of treatment-control pairs",
  main = "Hierarchy tier determining each pair",
  las = 1
)
graphics::legend(
  "topright",
  legend = rownames(tier_matrix),
  fill = c("#4C78A8", "#72B7B2", "#F2CF5B", "#B0B0B0"),
  inset = c(-0.36, 0),
  bty = "n"
)
grDevices::dev.off()

reproducibility_record <- c(
  paste("Analysis script:", normalizePath(script_path())),
  paste("Data file:", normalizePath(data_path)),
  paste("Data MD5:", unname(tools::md5sum(data_path))),
  paste("Dictionary file:", normalizePath(dictionary_path)),
  paste("Dictionary MD5:", unname(tools::md5sum(dictionary_path))),
  paste("Known randomization probability e:", known_e),
  "",
  capture.output(utils::sessionInfo())
)
writeLines(reproducibility_record, file.path(output_dir, "trial23_session_info.txt"))

report_lines <- c(
  "# Trial 23 fixed-horizon hierarchical win analysis",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Analysis specification",
  "",
  paste0("- Participants: N = ", N, "; rhTM n1 = ", n1, "; No rhTM n0 = ", n0, "."),
  paste0("- Known simple-randomization probability: e = ", known_e, "."),
  paste0("- Treatment coding: rhTM = 1; No rhTM = 0."),
  paste0("- YP_event_28d value denoting death/event: `", death_value, "`."),
  "- Hierarchy: 28-day survival > Day-3 change in SOFA > Day-3 change in D-dimer.",
  "- Lower Day-3 delta SOFA and lower Day-3 delta D-dimer are favorable.",
  "- Positive NB and WR/WO above 1 favor rhTM.",
  "- Primary inference uses the subject-level Hajek IF covariance estimator; no bootstrap is used.",
  "",
  "## Data validation",
  "",
  paste0("- Complete Day-28 event/time records: yes (", N, "/", N, ")."),
  paste0("- Day-28 survivors: ", sum(dat$death_event == 0L), "; all have Day-3 SOFA and D-dimer."),
  paste0("- Day-28 deaths: ", sum(dat$death_event == 1L), ". Missing Day-3 SOFA = ", sum(is.na(dat$YS_delta_SOFA_3d)), "; missing Day-3 D-dimer = ", sum(is.na(dat$YS_delta_Ddimer_3d)), "."),
  paste0("- Day-28 deaths by arm: rhTM = ", sum(dat$death_event == 1L & dat$A == 1L), "; No rhTM = ", sum(dat$death_event == 1L & dat$A == 0L), "."),
  "- No hierarchical outcome was imputed and no randomized subject was removed.",
  paste0("- Primary LADJ covariates: ", paste(primary_covariates, collapse = ", "), "."),
  paste0("- Optional X_TM sensitivity: pooled median = ", format_number(tm_median, 3), "; missing indicator count = ", sum(is.na(dat$X_TM)), "."),
  "",
  "## Pairwise hierarchy summary",
  "",
  markdown_table(pairwise_summary, 4),
  "",
  paste0(
    "Under the 2-point SOFA margin, ", reassigned_to_ddimer,
    " pairs that would have been decided by threshold-free SOFA were passed to D-dimer (",
    format_number(100 * reassigned_to_ddimer / (n1 * n0), 2), "% of all pairs)."
  ),
  "",
  "## Primary UNADJ and LADJ results",
  "",
  markdown_table(primary_results[c(
    "HCE_definition", "Method", "Measure", "Estimate_log_scale",
    "Estimate_report_scale", "SE", "CI_lower", "CI_upper", "P_value",
    "Estimated_variance", "Variance_reduction_percent", "Relative_efficiency"
  )], 4),
  "",
  "For WR and WO, SE and p-values are on the log scale; estimates and displayed confidence limits are exponentiated. For NB, inference is on the NB scale.",
  "",
  "## Primary efficiency comparison",
  "",
  markdown_table(primary_efficiency, 4),
  "",
  "Variance reduction is reported as computed and is not truncated at zero.",
  "",
  "## Baseline balance",
  "",
  markdown_table(baseline_balance, 4),
  "",
  "The standardized mean difference uses the pooled within-arm SD. Coefficients describe precision augmentation and are not causal effects.",
  "",
  "## Primary LADJ coefficients",
  "",
  markdown_table(primary_coefficients, 6),
  "",
  "## LADJ numerical diagnostics",
  "",
  markdown_table(primary_diagnostics, 6),
  "",
  "Covariate rescaling to pooled SDs produced the same fitted augmentation values up to the reported numerical tolerance.",
  "",
  "## Optional baseline-TM sensitivity",
  "",
  markdown_table(tm_sensitivity_results[c(
    "HCE_definition", "Method", "Measure", "Estimate_log_scale",
    "Estimate_report_scale", "SE", "CI_lower", "CI_upper", "P_value",
    "Variance_reduction_percent", "Relative_efficiency"
  )], 4),
  "",
  "This sensitivity analysis is separate from the primary result. X_TM was imputed using the pooled observed median without treatment or outcome information, and a missingness indicator was included.",
  "",
  "## Quality control",
  "",
  paste0("All ", nrow(quality_control), " prespecified quality-control checks passed."),
  "",
  markdown_table(quality_control[c("HCE_definition", "Check", "Observed", "Tolerance", "Passed")], 8),
  "",
  "## Interpretation cautions",
  "",
  "- This is a small trial; Wald inference relies on a first-order subject-level influence-function approximation.",
  "- The 2-point SOFA rule is a pairwise clinical relevance margin, not an individual responder definition.",
  "- Estimated variance reduction in this realized dataset need not be nonnegative even though LADJ has an asymptotic no-loss property under the stated projection conditions.",
  "- Clinical suitability also depends on whether the prespecified hierarchy and direction of the delta outcomes are defensible for the Trial 23 scientific question."
)
writeLines(report_lines, file.path(output_dir, "trial23_analysis_summary.md"))

cat("\nCompact primary analysis summary\n")
for (hce in pairwise_summary$HCE_definition) {
  ps <- pairwise_summary[pairwise_summary$HCE_definition == hce, ]
  rr <- primary_results[primary_results$HCE_definition == hce, ]
  ee <- primary_efficiency[primary_efficiency$HCE_definition == hce, ]
  cat("\n", hce, "\n", sep = "")
  cat("N / n1 / n0:", N, "/", n1, "/", n0, "\n")
  cat("Pair wins / losses / ties:", ps$Treatment_wins, "/", ps$Treatment_losses, "/", ps$Final_ties, "\n")
  cat("tau_w / tau_l / tau_tie:",
      paste(format_number(c(ps$tau_w, ps$tau_l, ps$tau_tie), 4), collapse = " / "), "\n")
  for (method in c("UNADJ", "LADJ")) {
    cat(method, ":\n", sep = "")
    method_rows <- rr[rr$Method == method, ]
    for (k in seq_len(nrow(method_rows))) {
      row <- method_rows[k, ]
      cat("  ", row$Measure, " = ", format_number(row$Estimate_report_scale, 4),
          "; SE = ", format_number(row$SE, 4),
          "; 95% CI [", format_number(row$CI_lower, 4), ", ", format_number(row$CI_upper, 4), "]",
          "; p = ", format_number(row$P_value, 4), "\n", sep = "")
    }
  }
  cat("Variance reduction NB / LWR / LWO (%): ",
      paste(format_number(ee$Variance_reduction_percent, 2), collapse = " / "), "\n", sep = "")
  cat("Relative efficiency NB / LWR / LWO: ",
      paste(format_number(ee$Relative_efficiency, 3), collapse = " / "), "\n", sep = "")
  cat("Tier counts: survival =", ps$Decided_by_survival,
      "; SOFA =", ps$Decided_by_SOFA,
      "; D-dimer =", ps$Decided_by_Ddimer,
      "; final ties =", ps$Final_ties, "\n")
  cat("Tier percentages: survival =", format_number(ps$Decided_by_survival_percent, 2),
      "%; SOFA =", format_number(ps$Decided_by_SOFA_percent, 2),
      "%; D-dimer =", format_number(ps$Decided_by_Ddimer_percent, 2),
      "%; final ties =", format_number(ps$Final_ties_percent, 2), "%\n", sep = "")
}
cat("\nAll outputs saved to:", normalizePath(output_dir), "\n")
