options(stringsAsFactors = FALSE, scipen = 999)
Sys.setenv(TZ = "America/New_York")

known_e <- 0.5
analysis_horizon <- "1 year"
primary_covariates <- c(
  "X_age_years",
  "X_sex_code",
  "X_hypertension",
  "X_dyslipidemia",
  "X_diabetes",
  "X_heart_failure",
  "X_prior_ami",
  "X_prior_pci",
  "X_prior_cabg"
)
sensitivity_label <- "LADJ sensitivity using all eligible complete baseline covariates"
outcome_variables <- c(
  "YP_death_1y",
  "YS_ami_1y",
  "YS_repeat_pci_1y",
  "YS_cabg_1y"
)

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", file_arg)))
  }
  normalizePath("analysis_trial88_win_ladj.R")
}

analysis_dir <- dirname(script_path())
data_path <- file.path(analysis_dir, "trial88.rds")
dictionary_path <- file.path(analysis_dir, "trial88_dictionary.csv")
output_dir <- file.path(analysis_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_analysis <- function(...) stop(..., call. = FALSE)
assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop_analysis(message)
}
write_csv <- function(x, filename) {
  utils::write.csv(x, file.path(output_dir, filename), row.names = FALSE, na = "")
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

assert_true(file.exists(data_path), paste("Missing RDS file:", data_path))
assert_true(file.exists(dictionary_path), paste("Missing dictionary:", dictionary_path))
dat <- as.data.frame(readRDS(data_path), check.names = FALSE)
dictionary <- utils::read.csv(dictionary_path, check.names = FALSE)
N <- nrow(dat)
assert_true(N == 200L, paste("Expected N=200 but found", N))
dat$participant_id <- seq_len(N)

required_variables <- c("Treatment", outcome_variables, primary_covariates)
missing_variables <- setdiff(required_variables, names(dat))
assert_true(
  length(missing_variables) == 0L,
  paste("Required variable(s) missing:", paste(missing_variables, collapse = ", "))
)
expected_treatment_levels <- c("Control", "Nutrition education")
treatment_levels <- unique(as.character(dat$Treatment))
assert_true(
  setequal(treatment_levels, expected_treatment_levels),
  paste("Unexpected treatment levels:", paste(treatment_levels, collapse = ", "))
)
dat$A <- as.integer(as.character(dat$Treatment) == "Nutrition education")
n1 <- sum(dat$A == 1L)
n0 <- sum(dat$A == 0L)
assert_true(n1 > 0L && n0 > 0L, "Both randomized arms must be present.")

for (variable in outcome_variables) {
  assert_true(!anyNA(dat[[variable]]), paste(variable, "contains missing values."))
  assert_true(
    identical(sort(unique(as.integer(dat[[variable]]))), 0:1),
    paste(variable, "must be binary with values 0 and 1.")
  )
  dictionary_row <- dictionary[dictionary$variable_name == variable, , drop = FALSE]
  assert_true(nrow(dictionary_row) == 1L, paste("Dictionary row missing for", variable))
  assert_true(
    grepl("1", dictionary_row$levels_or_range) &&
      grepl("Death|infarction|grafting|intervention", dictionary_row$short_explanation, ignore.case = TRUE),
    paste("Dictionary does not support event-occurrence coding for", variable)
  )
}
dat$revascularization_1y <- pmax(dat$YS_repeat_pci_1y, dat$YS_cabg_1y)
assert_true(
  all(dat$revascularization_1y %in% 0:1),
  "Derived repeat-revascularization status is not binary."
)
composite_from_components <- pmax(
  dat$YP_death_1y,
  dat$YS_ami_1y,
  dat$YS_repeat_pci_1y,
  dat$YS_cabg_1y
)
assert_true(
  "YS_composite_event_1y" %in% names(dat) &&
    identical(composite_from_components, dat$YS_composite_event_1y),
  "Adverse-event components do not reproduce YS_composite_event_1y."
)
assert_true(
  all(vapply(dat[primary_covariates], function(x) !anyNA(x), logical(1))),
  "A primary LADJ baseline covariate contains missing values."
)

dictionary_x <- dictionary[
  dictionary$variable_role == "Baseline covariate" &
    grepl("^X_", dictionary$variable_name) &
    dictionary$variable_name != "X_source_id",
  ,
  drop = FALSE
]
eligible_x <- dictionary_x$variable_name
missing_x_from_data <- setdiff(eligible_x, names(dat))
assert_true(
  length(missing_x_from_data) == 0L,
  paste("Eligible dictionary covariate(s) absent from data:", paste(missing_x_from_data, collapse = ", "))
)
assert_true(
  all(vapply(dat[eligible_x], function(x) !anyNA(x), logical(1))),
  "An eligible all-covariate sensitivity variable contains missing values."
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
    Used_in_primary = variable %in% required_variables,
    Eligible_for_sensitivity = variable %in% eligible_x,
    stringsAsFactors = FALSE
  )
}
data_audit <- do.call(rbind, lapply(audit_variables, audit_row))
dictionary_columns <- intersect(
  c("variable_name", "variable_role", "variable_type", "short_explanation"),
  names(dictionary)
)
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
write_csv(data_audit, "trial88_data_audit.csv")

cat("Trial 88 input audit\n")
cat("N =", N, "; n1 (Nutrition education) =", n1, "; n0 (Control) =", n0, "\n")
cat("Known randomization probability e =", known_e, "\n")
cat("Treatment levels:", paste(levels(dat$Treatment), collapse = ", "), "\n")
cat("All four adverse-event outcomes are complete binary 0/1 variables; dictionary coding defines 1 as occurrence.\n")
cat("Primary LADJ covariates:", paste(primary_covariates, collapse = ", "), "\n")
cat("Eligible complete baseline X variables for sensitivity:", paste(eligible_x, collapse = ", "), "\n\n")

make_sensitivity_features <- function(data, eligible_names) {
  procedure_names <- c(
    "X_balloon_angioplasty",
    "X_stent",
    "X_pharmacologic_stent"
  )
  assert_true(
    all(procedure_names %in% eligible_names),
    "Expected baseline procedure indicators are not all eligible."
  )
  procedure_sum <- rowSums(data[procedure_names])
  assert_true(
    all(procedure_sum == 1L),
    "Procedure indicators are not mutually exclusive and exhaustive."
  )
  procedure_level <- apply(data[procedure_names], 1L, function(row) procedure_names[which(row == 1L)])
  procedure_factor <- factor(procedure_level, levels = procedure_names)
  procedure_dummies <- stats::model.matrix(~ procedure_factor)[, -1L, drop = FALSE]
  colnames(procedure_dummies) <- c(
    "X_procedure_stent_vs_balloon",
    "X_procedure_pharmacologic_stent_vs_balloon"
  )

  smoking_levels <- sort(unique(data$X_smoking_baseline))
  assert_true(identical(smoking_levels, 0:2), "Unexpected X_smoking_baseline levels.")
  smoking_factor <- factor(data$X_smoking_baseline, levels = smoking_levels)
  smoking_dummies <- stats::model.matrix(~ smoking_factor)[, -1L, drop = FALSE]
  colnames(smoking_dummies) <- paste0("X_smoking_baseline_", smoking_levels[-1L], "_vs_0")

  retained_raw <- setdiff(
    eligible_names,
    c(procedure_names, "X_smoking_baseline")
  )
  X <- cbind(
    as.data.frame(data[retained_raw]),
    procedure_dummies,
    smoking_dummies
  )
  list(
    X = as.data.frame(X, check.names = FALSE),
    retained_raw = retained_raw,
    procedure_names = procedure_names,
    procedure_reference = "X_balloon_angioplasty",
    smoking_reference = "0"
  )
}

sensitivity_features <- make_sensitivity_features(dat, eligible_x)
primary_X <- dat[primary_covariates]
sensitivity_X <- sensitivity_features$X

raw_smoking_factor <- factor(dat$X_smoking_baseline, levels = 0:2)
raw_smoking_dummies <- stats::model.matrix(~ raw_smoking_factor)[, -1L, drop = FALSE]
colnames(raw_smoking_dummies) <- c("X_smoking_baseline_1_vs_0", "X_smoking_baseline_2_vs_0")
raw_all_candidate <- cbind(
  dat[setdiff(eligible_x, "X_smoking_baseline")],
  raw_smoking_dummies
)
raw_all_centered <- scale(as.matrix(raw_all_candidate), center = TRUE, scale = FALSE)
raw_all_gram <- crossprod(raw_all_centered * (dat$A - known_e))
raw_all_eigen <- eigen(raw_all_gram, symmetric = TRUE, only.values = TRUE)$values
sensitivity_centered <- scale(as.matrix(sensitivity_X), center = TRUE, scale = FALSE)
max_sensitivity_span_error <- max(abs(qr.resid(qr(sensitivity_centered), raw_all_centered)))
raw_all_feature_audit <- data.frame(
  Candidate_coding = "All dictionary X columns; smoking dummy-coded; procedure indicators unreparameterized",
  Number_columns = ncol(raw_all_candidate),
  Gram_rank = qr(raw_all_gram)$rank,
  Smallest_eigenvalue = min(raw_all_eigen),
  Largest_eigenvalue = max(raw_all_eigen),
  Condition_number = if (min(raw_all_eigen) > 0) max(raw_all_eigen) / min(raw_all_eigen) else Inf,
  Exact_relation = paste(
    "X_balloon_angioplasty + X_stent + X_pharmacologic_stent = 1;",
    "centered columns sum to zero"
  ),
  Resolution = paste(
    "Reference-code the three-level procedure category with balloon angioplasty as reference;",
    "this preserves its full two-dimensional centered span"
  ),
  Max_raw_column_residual_from_final_coded_span = max_sensitivity_span_error,
  stringsAsFactors = FALSE
)
write_csv(raw_all_feature_audit, "trial88_sensitivity_feature_audit.csv")

make_pairwise_comparisons <- function(data) {
  treated_rows <- which(data$A == 1L)
  control_rows <- which(data$A == 0L)
  n_t <- length(treated_rows)
  n_c <- length(control_rows)
  win <- matrix(0L, nrow = n_t, ncol = n_c)
  loss <- matrix(0L, nrow = n_t, ncol = n_c)
  tie <- matrix(0L, nrow = n_t, ncol = n_c)
  decision_tier <- matrix(NA_character_, nrow = n_t, ncol = n_c)

  tr_death <- matrix(data$YP_death_1y[treated_rows], nrow = n_t, ncol = n_c)
  co_death <- matrix(data$YP_death_1y[control_rows], nrow = n_t, ncol = n_c, byrow = TRUE)
  death_win <- tr_death == 0L & co_death == 1L
  death_loss <- tr_death == 1L & co_death == 0L
  both_die <- tr_death == 1L & co_death == 1L
  both_survive <- tr_death == 0L & co_death == 0L
  win[death_win] <- 1L
  loss[death_loss] <- 1L
  tie[both_die] <- 1L
  decision_tier[death_win | death_loss] <- "1-year survival"
  decision_tier[both_die] <- "Final tie at death"

  tr_ami <- matrix(data$YS_ami_1y[treated_rows], nrow = n_t, ncol = n_c)
  co_ami <- matrix(data$YS_ami_1y[control_rows], nrow = n_t, ncol = n_c, byrow = TRUE)
  ami_win <- both_survive & tr_ami == 0L & co_ami == 1L
  ami_loss <- both_survive & tr_ami == 1L & co_ami == 0L
  pass_ami <- both_survive & tr_ami == co_ami
  win[ami_win] <- 1L
  loss[ami_loss] <- 1L
  decision_tier[ami_win | ami_loss] <- "AMI-free status"

  tr_revasc <- matrix(data$revascularization_1y[treated_rows], nrow = n_t, ncol = n_c)
  co_revasc <- matrix(data$revascularization_1y[control_rows], nrow = n_t, ncol = n_c, byrow = TRUE)
  revasc_win <- pass_ami & tr_revasc == 0L & co_revasc == 1L
  revasc_loss <- pass_ami & tr_revasc == 1L & co_revasc == 0L
  final_tie <- pass_ami & tr_revasc == co_revasc
  win[revasc_win] <- 1L
  loss[revasc_loss] <- 1L
  tie[final_tie] <- 1L
  decision_tier[revasc_win | revasc_loss] <- "Repeat-revascularization-free status"
  decision_tier[final_tie] <- "Final tie at revascularization"

  assert_true(all(win + loss + tie == 1L), "Pairwise win/loss/tie partition failed.")
  assert_true(!anyNA(decision_tier), "At least one pair lacks a decision tier.")
  pair_index <- expand.grid(
    treated_row = treated_rows,
    control_row = control_rows,
    KEEP.OUT.ATTRS = FALSE
  )
  pairs <- data.frame(
    treated_id = data$participant_id[pair_index$treated_row],
    control_id = data$participant_id[pair_index$control_row],
    win_ij = as.integer(win),
    loss_ij = as.integer(loss),
    tie_ij = as.integer(tie),
    decision_tier = as.character(decision_tier),
    stringsAsFactors = FALSE
  )
  list(
    win = win,
    loss = loss,
    tie = tie,
    decision_tier = decision_tier,
    pairs = pairs,
    treated_rows = treated_rows,
    control_rows = control_rows
  )
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
  assert_true(tau_w > 0 && tau_l > 0, "Win/loss probability must be positive for LWR.")
  delta <- tau_w - tau_l
  assert_true(abs(delta) < 1, "Absolute net benefit must be less than one for LWO.")
  c(
    NB = delta,
    WR = log(tau_w) - log(tau_l),
    WO = log1p(delta) - log1p(-delta)
  )
}

make_inference_rows <- function(method, tau, psi_components) {
  values <- measure_values(tau[["w"]], tau[["l"]])
  psi_components <- as.matrix(psi_components)
  centered <- sweep(psi_components, 2L, colMeans(psi_components), "-")
  sigma <- crossprod(centered) / N
  rows <- lapply(names(values), function(measure) {
    gradient <- measure_gradient(measure, tau[["w"]], tau[["l"]])
    variance <- as.numeric(crossprod(gradient, sigma %*% gradient) / N)
    assert_true(is.finite(variance) && variance >= -1e-14, "Invalid estimated variance.")
    variance <- max(variance, 0)
    se <- sqrt(variance)
    estimate_log <- unname(values[[measure]])
    ci_log <- estimate_log + stats::qnorm(c(0.025, 0.975)) * se
    p_value <- 2 * stats::pnorm(-abs(estimate_log / se))
    if (measure == "NB") {
      report_estimate <- estimate_log
      report_ci <- ci_log
      log_output <- NA_real_
    } else {
      report_estimate <- exp(estimate_log)
      report_ci <- exp(ci_log)
      log_output <- estimate_log
    }
    data.frame(
      Method = method,
      Measure = measure,
      Estimate_log_scale = log_output,
      Estimate_report_scale = report_estimate,
      SE = se,
      CI_lower = report_ci[[1L]],
      CI_upper = report_ci[[2L]],
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

fit_ladj <- function(scores, X, covariate_names, method_label) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  assert_true(!anyNA(X), paste("Missing LADJ feature in", method_label))
  assert_true(ncol(X) == length(covariate_names), "Feature names do not match LADJ matrix.")
  X_center <- sweep(X, 2L, colMeans(X), "-")
  colnames(X_center) <- covariate_names
  Z <- dat$A - known_e
  design <- X_center * Z
  gram <- crossprod(design)
  rank <- qr(gram)$rank
  eigenvalues <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
  smallest_eigenvalue <- min(eigenvalues)
  largest_eigenvalue <- max(eigenvalues)
  assert_true(rank == ncol(X_center), paste("Singular LADJ Gram matrix in", method_label))
  assert_true(
    smallest_eigenvalue > .Machine$double.eps * largest_eigenvalue * ncol(X_center),
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
  scaled_eigenvalues <- eigen(gram_scaled, symmetric = TRUE, only.values = TRUE)$values
  beta_scaled <- solve(gram_scaled, crossprod(design_scaled, scores$psi))
  fitted_scaled <- X_scaled %*% beta_scaled
  covariance <- crossprod(sweep(psi_adjusted, 2L, colMeans(psi_adjusted), "-")) / N

  list(
    tau = tau_adjusted,
    psi = psi_adjusted,
    covariance = covariance,
    coefficients = data.frame(
      Method = method_label,
      Covariate = covariate_names,
      beta_w = beta[, "beta_w"],
      beta_l = beta[, "beta_l"],
      stringsAsFactors = FALSE
    ),
    diagnostics = data.frame(
      Method = method_label,
      Number_features = ncol(X_center),
      Gram_rank = rank,
      Smallest_eigenvalue = smallest_eigenvalue,
      Largest_eigenvalue = largest_eigenvalue,
      Condition_number = largest_eigenvalue / smallest_eigenvalue,
      Scaled_condition_number = max(scaled_eigenvalues) / min(scaled_eigenvalues),
      Max_centering_identity_error = max(abs(identity_left - identity_right)),
      Max_tau_formula_error = max(abs(direct_correction - imbalance_correction)),
      Max_scaled_fitted_value_error = max(abs(fitted_b - fitted_scaled)),
      Min_adjusted_covariance_eigenvalue = min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values),
      stringsAsFactors = FALSE
    )
  )
}

pooled_within_arm_sd <- function(x, A) {
  n_t <- sum(A == 1L)
  n_c <- sum(A == 0L)
  sqrt(((n_t - 1) * stats::var(x[A == 1L]) + (n_c - 1) * stats::var(x[A == 0L])) /
         (n_t + n_c - 2))
}
make_balance <- function(X, names) {
  X <- as.data.frame(X, check.names = FALSE)
  do.call(rbind, lapply(names, function(variable) {
    x <- X[[variable]]
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
}

pair_object <- make_pairwise_comparisons(dat)
scores <- construct_scores(pair_object)
primary_fit <- fit_ladj(scores, primary_X, primary_covariates, "LADJ")
sensitivity_fit <- fit_ladj(
  scores,
  sensitivity_X,
  colnames(sensitivity_X),
  sensitivity_label
)
unadjusted_inference <- make_inference_rows("UNADJ", scores$tau, scores$psi)
ladj_inference <- make_inference_rows("LADJ", primary_fit$tau, primary_fit$psi)
sensitivity_inference <- make_inference_rows(
  sensitivity_label,
  sensitivity_fit$tau,
  sensitivity_fit$psi
)
results <- rbind(
  unadjusted_inference$rows,
  ladj_inference$rows,
  sensitivity_inference$rows
)
unadjusted_variance <- setNames(
  unadjusted_inference$rows$Estimated_variance,
  unadjusted_inference$rows$Measure
)
make_efficiency <- function(method_rows) {
  adjusted_variance <- setNames(method_rows$Estimated_variance, method_rows$Measure)
  data.frame(
    Method = unique(method_rows$Method),
    Measure = names(unadjusted_variance),
    Unadjusted_variance = as.numeric(unadjusted_variance),
    Adjusted_variance = as.numeric(adjusted_variance[names(unadjusted_variance)]),
    Variance_reduction_percent = 100 * (1 - adjusted_variance[names(unadjusted_variance)] /
                                          unadjusted_variance),
    Relative_efficiency = unadjusted_variance /
      adjusted_variance[names(unadjusted_variance)],
    stringsAsFactors = FALSE
  )
}
primary_efficiency <- make_efficiency(ladj_inference$rows)
sensitivity_efficiency <- make_efficiency(sensitivity_inference$rows)
results$Variance_reduction_percent <- 0
results$Relative_efficiency <- 1
for (eff in list(primary_efficiency, sensitivity_efficiency)) {
  for (row_id in seq_len(nrow(eff))) {
    idx <- results$Method == eff$Method[[row_id]] & results$Measure == eff$Measure[[row_id]]
    results$Variance_reduction_percent[idx] <- eff$Variance_reduction_percent[[row_id]]
    results$Relative_efficiency[idx] <- eff$Relative_efficiency[[row_id]]
  }
}

pairs <- pair_object$pairs
N_pairs <- nrow(pairs)
pairwise_summary <- data.frame(
  N = N,
  n1 = n1,
  n0 = n0,
  N_pairs = N_pairs,
  Treatment_wins = sum(pairs$win_ij),
  Treatment_losses = sum(pairs$loss_ij),
  Final_ties = sum(pairs$tie_ij),
  tau_w = mean(pairs$win_ij),
  tau_l = mean(pairs$loss_ij),
  tau_tie = mean(pairs$tie_ij),
  Decided_by_death = sum(pairs$decision_tier == "1-year survival"),
  Decided_by_death_percent = 100 * mean(pairs$decision_tier == "1-year survival"),
  Decided_by_AMI = sum(pairs$decision_tier == "AMI-free status"),
  Decided_by_AMI_percent = 100 * mean(pairs$decision_tier == "AMI-free status"),
  Decided_by_revascularization = sum(pairs$decision_tier == "Repeat-revascularization-free status"),
  Decided_by_revascularization_percent = 100 * mean(pairs$decision_tier == "Repeat-revascularization-free status"),
  Final_ties_percent = 100 * mean(pairs$tie_ij),
  Final_ties_at_death = sum(pairs$decision_tier == "Final tie at death"),
  Final_ties_at_revascularization = sum(pairs$decision_tier == "Final tie at revascularization"),
  stringsAsFactors = FALSE
)

primary_balance <- make_balance(primary_X, primary_covariates)
sensitivity_balance <- make_balance(sensitivity_X, colnames(sensitivity_X))
subject_scores <- data.frame(
  participant_id = dat$participant_id,
  Treatment = as.character(dat$Treatment),
  A = dat$A,
  psi_w = scores$psi[, "w"],
  psi_l = scores$psi[, "l"],
  psi_w_LADJ = primary_fit$psi[, "w"],
  psi_l_LADJ = primary_fit$psi[, "l"],
  stringsAsFactors = FALSE
)

projection_reduction <- function(unadjusted_score, adjusted_score) {
  100 * (1 - stats::var(adjusted_score) / stats::var(unadjusted_score))
}
projection_strength <- rbind(
  data.frame(
    Score = c("psi_w", "psi_l", "psi_w_plus_psi_l", "psi_w_minus_psi_l"),
    Interpretation = c(
      "Elementary win IF",
      "Elementary loss IF",
      "Non-tie/activity direction",
      "NB direction"
    ),
    Fixed_gradient_projection_reduction_percent = c(
      projection_reduction(scores$psi[, "w"], primary_fit$psi[, "w"]),
      projection_reduction(scores$psi[, "l"], primary_fit$psi[, "l"]),
      projection_reduction(rowSums(scores$psi), rowSums(primary_fit$psi)),
      projection_reduction(
        scores$psi[, "w"] - scores$psi[, "l"],
        primary_fit$psi[, "w"] - primary_fit$psi[, "l"]
      )
    ),
    stringsAsFactors = FALSE
  ),
  do.call(rbind, lapply(c("NB", "WR", "WO"), function(measure) {
    gradient <- measure_gradient(measure, scores$tau[["w"]], scores$tau[["l"]])
    data.frame(
      Score = paste0(measure, "_IF_at_unadjusted_gradient"),
      Interpretation = paste(measure, "score using one common fixed gradient"),
      Fixed_gradient_projection_reduction_percent = projection_reduction(
        as.vector(scores$psi %*% gradient),
        as.vector(primary_fit$psi %*% gradient)
      ),
      stringsAsFactors = FALSE
    )
  }))
)
score_matrix <- cbind(
  psi_w = scores$psi[, "w"],
  psi_l = scores$psi[, "l"],
  NB_IF = scores$psi[, "w"] - scores$psi[, "l"]
)
lwr_gradient_unadjusted <- measure_gradient("WR", scores$tau[["w"]], scores$tau[["l"]])
score_matrix <- cbind(
  score_matrix,
  LWR_IF = as.vector(scores$psi %*% lwr_gradient_unadjusted)
)
score_correlations_matrix <- stats::cor(as.matrix(primary_X), score_matrix)
score_correlations <- do.call(rbind, lapply(seq_len(nrow(score_correlations_matrix)), function(i) {
  data.frame(
    Covariate = rownames(score_correlations_matrix)[[i]],
    Score = colnames(score_correlations_matrix),
    Correlation = as.numeric(score_correlations_matrix[i, ]),
    stringsAsFactors = FALSE
  )
}))

add_qc <- function(check, observed, tolerance, passed, details = "") {
  data.frame(
    Check = check,
    Observed = observed,
    Tolerance = tolerance,
    Passed = passed,
    Details = details,
    stringsAsFactors = FALSE
  )
}
placement_means <- c(
  mean(scores$placement$g_w_T),
  mean(scores$placement$g_l_T),
  mean(scores$placement$g_w_C),
  mean(scores$placement$g_l_C)
)
primary_covariance_error <- max(abs(primary_fit$covariance - t(primary_fit$covariance)))
sensitivity_covariance_error <- max(abs(sensitivity_fit$covariance - t(sensitivity_fit$covariance)))
quality_control <- rbind(
  add_qc("Number of pairs equals n1*n0", N_pairs, n1 * n0, N_pairs == n1 * n0),
  add_qc("Each pair has win+loss+tie=1", max(abs(pairs$win_ij + pairs$loss_ij + pairs$tie_ij - 1)), 0, all(pairs$win_ij + pairs$loss_ij + pairs$tie_ij == 1)),
  add_qc("Pair probabilities sum to one", abs(sum(scores$tau) + mean(pairs$tie_ij) - 1), 1e-12, abs(sum(scores$tau) + mean(pairs$tie_ij) - 1) < 1e-12),
  add_qc("Placement means are zero", max(abs(placement_means)), 1e-12, max(abs(placement_means)) < 1e-12),
  add_qc("Mean elementary IFs are zero", max(abs(colMeans(scores$psi))), 1e-12, max(abs(colMeans(scores$psi))) < 1e-12),
  add_qc("Known e is exactly 0.5", known_e, 0.5, identical(known_e, 0.5)),
  add_qc("No intercept in LADJ", 0, 0, TRUE),
  add_qc("Primary LADJ uses exact prespecified covariates", length(primary_covariates), 9, identical(colnames(primary_X), primary_covariates)),
  add_qc("Primary adjusted tau formulas agree", primary_fit$diagnostics$Max_tau_formula_error, 1e-12, primary_fit$diagnostics$Max_tau_formula_error < 1e-12),
  add_qc("Sensitivity adjusted tau formulas agree", sensitivity_fit$diagnostics$Max_tau_formula_error, 1e-12, sensitivity_fit$diagnostics$Max_tau_formula_error < 1e-12),
  add_qc("Primary Gram matrix is full rank", primary_fit$diagnostics$Gram_rank, ncol(primary_X), primary_fit$diagnostics$Gram_rank == ncol(primary_X)),
  add_qc("Sensitivity final coding is full rank", sensitivity_fit$diagnostics$Gram_rank, ncol(sensitivity_X), sensitivity_fit$diagnostics$Gram_rank == ncol(sensitivity_X)),
  add_qc("Raw all-X rank deficiency was detected", qr(raw_all_gram)$rank, ncol(raw_all_candidate) - 1L, qr(raw_all_gram)$rank == ncol(raw_all_candidate) - 1L, raw_all_feature_audit$Exact_relation),
  add_qc("Final sensitivity coding preserves raw all-X centered span", max_sensitivity_span_error, 1e-12, max_sensitivity_span_error < 1e-12),
  add_qc("Primary covariance is symmetric", primary_covariance_error, 1e-12, primary_covariance_error < 1e-12),
  add_qc("Sensitivity covariance is symmetric", sensitivity_covariance_error, 1e-12, sensitivity_covariance_error < 1e-12),
  add_qc("Primary covariance is positive semidefinite", primary_fit$diagnostics$Min_adjusted_covariance_eigenvalue, -1e-12, primary_fit$diagnostics$Min_adjusted_covariance_eigenvalue >= -1e-12),
  add_qc("Sensitivity covariance is positive semidefinite", sensitivity_fit$diagnostics$Min_adjusted_covariance_eigenvalue, -1e-12, sensitivity_fit$diagnostics$Min_adjusted_covariance_eigenvalue >= -1e-12),
  add_qc("All four outcomes are complete", sum(is.na(dat[outcome_variables])), 0, !anyNA(dat[outcome_variables])),
  add_qc("All four outcomes use 0/1 adverse-event coding", 1, 1, all(vapply(dat[outcome_variables], function(x) identical(sort(unique(x)), 0:1), logical(1)))),
  add_qc("Adverse components reproduce composite event", sum(composite_from_components != dat$YS_composite_event_1y), 0, identical(composite_from_components, dat$YS_composite_event_1y)),
  add_qc("All randomized subjects retained", length(unique(c(pairs$treated_id, pairs$control_id))), N, length(unique(c(pairs$treated_id, pairs$control_id))) == N),
  add_qc("tau_w is positive for LWR", scores$tau[["w"]], 0, scores$tau[["w"]] > 0),
  add_qc("tau_l is positive for LWR", scores$tau[["l"]], 0, scores$tau[["l"]] > 0)
)
assert_true(all(quality_control$Passed), "At least one quality-control check failed.")

primary_results <- results[results$Method %in% c("UNADJ", "LADJ"), ]
sensitivity_results <- results[results$Method == sensitivity_label, ]
write_csv(pairs, "trial88_pairwise_comparisons.csv")
write_csv(pairwise_summary, "trial88_pairwise_summary.csv")
write_csv(primary_results, "trial88_unadjusted_ladj_results.csv")
write_csv(primary_efficiency, "trial88_efficiency_comparison.csv")
write_csv(primary_balance, "trial88_baseline_balance.csv")
write_csv(primary_fit$coefficients, "trial88_ladj_coefficients.csv")
write_csv(primary_fit$diagnostics, "trial88_ladj_diagnostics.csv")
write_csv(subject_scores, "trial88_subject_scores.csv")
write_csv(projection_strength, "trial88_projection_strength_diagnostics.csv")
write_csv(score_correlations, "trial88_covariate_score_correlations.csv")
write_csv(quality_control, "trial88_quality_control.csv")
write_csv(sensitivity_results, "trial88_all_covariate_sensitivity_results.csv")
write_csv(sensitivity_efficiency, "trial88_all_covariate_sensitivity_efficiency.csv")
write_csv(sensitivity_balance, "trial88_all_covariate_balance.csv")
write_csv(sensitivity_fit$coefficients, "trial88_all_covariate_sensitivity_coefficients.csv")
write_csv(sensitivity_fit$diagnostics, "trial88_all_covariate_sensitivity_diagnostics.csv")

min_tau <- min(scores$tau)
lwr_stability <- if (min_tau >= 0.05) {
  "Both elementary probabilities exceed 0.05; LWR is defined, although sparse events still warrant caution."
} else {
  "At least one elementary probability is below 0.05; LWR is mathematically defined but potentially unstable."
}
sensitivity_condition_flag <- if (sensitivity_fit$diagnostics$Scaled_condition_number < 100) {
  "The final SD-scaled sensitivity Gram matrix is acceptably conditioned."
} else {
  "The final SD-scaled sensitivity Gram matrix remains poorly conditioned; interpret sensitivity cautiously."
}

open_plot_device(
  file.path(output_dir, "trial88_primary_effect_estimates.png"),
  width = 1750,
  height = 750,
  res = 170
)
graphics::par(mfrow = c(1, 3), mar = c(5, 7, 4, 2), oma = c(0, 0, 2, 0))
method_colors <- c(UNADJ = "#5B6770", LADJ = "#1F77B4")
for (measure in c("NB", "WR", "WO")) {
  plot_data <- primary_results[primary_results$Measure == measure, ]
  y <- rev(seq_len(nrow(plot_data)))
  xlim <- range(c(plot_data$CI_lower, plot_data$CI_upper))
  padding <- diff(xlim) * 0.1
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
  graphics::segments(plot_data$CI_lower, y, plot_data$CI_upper, y, lwd = 2, col = method_colors[plot_data$Method])
  graphics::points(plot_data$Estimate_report_scale, y, pch = 19, cex = 1.1, col = method_colors[plot_data$Method])
  graphics::axis(2, at = y, labels = plot_data$Method, las = 1)
}
graphics::mtext("Trial 88 primary estimates with 95% Wald confidence intervals", outer = TRUE, cex = 1.15)
grDevices::dev.off()

tier_percentages <- c(
  "1-year survival" = pairwise_summary$Decided_by_death_percent,
  "AMI-free status" = pairwise_summary$Decided_by_AMI_percent,
  "Revascularization-free status" = pairwise_summary$Decided_by_revascularization_percent,
  "Final tie" = pairwise_summary$Final_ties_percent
)
open_plot_device(
  file.path(output_dir, "trial88_pair_decision_tiers.png"),
  width = 1300,
  height = 850,
  res = 170
)
graphics::par(mar = c(5, 5, 4, 12), xpd = TRUE)
graphics::barplot(
  matrix(tier_percentages, ncol = 1L),
  names.arg = "Trial 88 hierarchy",
  col = c("#4C78A8", "#72B7B2", "#F2CF5B", "#B0B0B0"),
  ylim = c(0, 100),
  ylab = "Percentage of treatment-control pairs",
  main = "Hierarchy tier determining each pair",
  las = 1
)
graphics::legend(
  "topright",
  legend = names(tier_percentages),
  fill = c("#4C78A8", "#72B7B2", "#F2CF5B", "#B0B0B0"),
  inset = c(-0.48, 0),
  bty = "n"
)
grDevices::dev.off()

event_counts <- do.call(rbind, lapply(outcome_variables, function(variable) {
  data.frame(
    Outcome = variable,
    Nutrition_education_events = sum(dat[[variable]][dat$A == 1L]),
    Control_events = sum(dat[[variable]][dat$A == 0L]),
    Total_events = sum(dat[[variable]]),
    stringsAsFactors = FALSE
  )
}))
write_csv(event_counts, "trial88_event_counts.csv")

report_lines <- c(
  "# Trial 88 fixed-horizon hierarchical win analysis",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Analysis specification",
  "",
  paste0("- Participants: N = ", N, "; Nutrition education n1 = ", n1, "; Control n0 = ", n0, "."),
  paste0("- Known simple-randomization probability: e = ", known_e, "; observed n1/N was not substituted."),
  "- Treatment coding: Nutrition education = 1; Control = 0.",
  paste0("- Fixed horizon: ", analysis_horizon, "."),
  "- Hierarchy: 1-year survival > AMI-free status > repeat-revascularization-free status.",
  "- Two deaths are tied overall because death times are unavailable.",
  "- Positive NB and WR/WO above 1 favor Nutrition education.",
  "- Primary inference uses the subject-level Hajek IF covariance estimator; no bootstrap is used.",
  "",
  "## Data validation and event sparsity",
  "",
  "All four outcomes are complete binary variables, and the dictionary defines 1 as occurrence of the adverse event.",
  "",
  markdown_table(event_counts, 0),
  "",
  "All nine primary covariates are baseline, pretreatment variables in the dictionary and have zero missingness.",
  "",
  "## Pairwise hierarchy summary",
  "",
  markdown_table(pairwise_summary, 4),
  "",
  paste0("Final tie probability = ", format_number(pairwise_summary$tau_tie, 4), ". ", lwr_stability),
  "",
  "## Primary UNADJ and LADJ results",
  "",
  markdown_table(primary_results[c(
    "Method", "Measure", "Estimate_log_scale", "Estimate_report_scale",
    "SE", "CI_lower", "CI_upper", "P_value", "Estimated_variance",
    "Variance_reduction_percent", "Relative_efficiency", "tau_w", "tau_l", "tau_tie"
  )], 5),
  "",
  "For WR and WO, SE and p-values are on the log scale; reported estimates and confidence limits are exponentiated. Final adjusted measures are exact transformations of adjusted tau_w and tau_l.",
  "",
  "## Primary efficiency comparison",
  "",
  markdown_table(primary_efficiency, 5),
  "",
  "Variance reduction is reported as computed and is not truncated at zero.",
  "",
  "## Primary linear projection strength",
  "",
  markdown_table(projection_strength, 4),
  "",
  paste(
    "These fixed-gradient reductions isolate how much of each IF score is removed by the linear projection.",
    "They differ in finite samples from comparing each estimator's final plug-in variance because WR gradients are evaluated at different unadjusted and adjusted elementary estimates."
  ),
  "",
  "## Primary baseline balance",
  "",
  markdown_table(primary_balance, 4),
  "",
  "The standardized mean difference uses the pooled within-arm SD.",
  "",
  "## Primary LADJ coefficients",
  "",
  markdown_table(primary_fit$coefficients, 7),
  "",
  "## Primary LADJ diagnostics",
  "",
  markdown_table(primary_fit$diagnostics, 7),
  "",
  "## All-covariate sensitivity coding",
  "",
  markdown_table(raw_all_feature_audit, 5),
  "",
  paste0(
    "Smoking was dummy-coded with 0 as reference. The three mutually exclusive procedure indicators were treated as one categorical variable; ",
    sensitivity_features$procedure_reference,
    " was the reference. This removes only a redundant centered column and preserves all procedure-category information."
  ),
  "",
  "## All-covariate sensitivity results",
  "",
  markdown_table(sensitivity_results[c(
    "Method", "Measure", "Estimate_log_scale", "Estimate_report_scale",
    "SE", "CI_lower", "CI_upper", "P_value", "Variance_reduction_percent",
    "Relative_efficiency", "tau_w", "tau_l", "tau_tie"
  )], 5),
  "",
  markdown_table(sensitivity_fit$diagnostics, 7),
  "",
  sensitivity_condition_flag,
  "",
  "## Quality control",
  "",
  paste0("All ", nrow(quality_control), " prespecified quality-control checks passed."),
  "",
  markdown_table(quality_control[c("Check", "Observed", "Tolerance", "Passed", "Details")], 8),
  "",
  "## Interpretation cautions",
  "",
  "- Cardiovascular events are sparse and most treatment-control pairs are expected to remain tied; the effective information is determined by subjects, not by treating all pairs as independent.",
  "- Wald inference relies on the first-order subject-level influence-function approximation.",
  "- The all-covariate analysis is a sensitivity analysis, not the primary analysis; it uses many predictors relative to N=200.",
  "- The observed variance reduction may be negative even though the projection no-loss statement is asymptotic."
)
writeLines(report_lines, file.path(output_dir, "trial88_analysis_summary.md"))

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
writeLines(reproducibility_record, file.path(output_dir, "trial88_session_info.txt"))

cat("\nCompact Trial 88 primary analysis summary\n")
cat("N / n1 / n0:", N, "/", n1, "/", n0, "\n")
cat("Pairs:", N_pairs, "; wins / losses / ties:",
    pairwise_summary$Treatment_wins, "/", pairwise_summary$Treatment_losses, "/", pairwise_summary$Final_ties, "\n")
cat("tau_w / tau_l / tau_tie:",
    paste(format_number(c(pairwise_summary$tau_w, pairwise_summary$tau_l, pairwise_summary$tau_tie), 5), collapse = " / "), "\n")
cat("Tier counts: death =", pairwise_summary$Decided_by_death,
    "; AMI =", pairwise_summary$Decided_by_AMI,
    "; revascularization =", pairwise_summary$Decided_by_revascularization,
    "; final ties =", pairwise_summary$Final_ties, "\n")
cat("Tier percentages: death =", format_number(pairwise_summary$Decided_by_death_percent, 2),
    "%; AMI =", format_number(pairwise_summary$Decided_by_AMI_percent, 2),
    "%; revascularization =", format_number(pairwise_summary$Decided_by_revascularization_percent, 2),
    "%; final ties =", format_number(pairwise_summary$Final_ties_percent, 2), "%\n", sep = "")
for (method in c("UNADJ", "LADJ")) {
  cat(method, ":\n", sep = "")
  method_rows <- primary_results[primary_results$Method == method, ]
  for (k in seq_len(nrow(method_rows))) {
    row <- method_rows[k, ]
    cat("  ", row$Measure, " = ", format_number(row$Estimate_report_scale, 5),
        "; SE = ", format_number(row$SE, 5),
        "; 95% CI [", format_number(row$CI_lower, 5), ", ", format_number(row$CI_upper, 5), "]",
        "; p = ", format_number(row$P_value, 5), "\n", sep = "")
  }
}
cat("Adjusted tau_w / tau_l:",
    paste(format_number(primary_fit$tau, 5), collapse = " / "), "\n")
cat("Variance reduction NB / LWR / LWO (%):",
    paste(format_number(primary_efficiency$Variance_reduction_percent, 2), collapse = " / "), "\n")
cat("Relative efficiency NB / LWR / LWO:",
    paste(format_number(primary_efficiency$Relative_efficiency, 3), collapse = " / "), "\n")
cat("Primary Gram condition number (raw / SD-scaled):",
    format_number(primary_fit$diagnostics$Condition_number, 2), "/",
    format_number(primary_fit$diagnostics$Scaled_condition_number, 2), "\n")
cat("Sensitivity Gram condition number (raw / SD-scaled):",
    format_number(sensitivity_fit$diagnostics$Condition_number, 2), "/",
    format_number(sensitivity_fit$diagnostics$Scaled_condition_number, 2), "\n")
cat("LWR stability:", lwr_stability, "\n")
cat("All outputs saved to:", normalizePath(output_dir), "\n")
