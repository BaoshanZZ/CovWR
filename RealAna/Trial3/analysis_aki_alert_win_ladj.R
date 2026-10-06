options(stringsAsFactors = FALSE, scipen = 999)
Sys.setenv(TZ = "America/New_York")

known_e <- 0.5
analysis_horizon_days <- 14L

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", file_arg)))
  }
  normalizePath("RealAna/Trial3/analysis_aki_alert_win_ladj.R")
}

project_dir <- dirname(dirname(dirname(script_path())))

stop_analysis <- function(...) stop(..., call. = FALSE)
assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop_analysis(message)
}
write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
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

trial_configs <- list(
  Trial3 = list(
    trial_id = 3L,
    data_file = file.path(project_dir, "RealAna", "Trial3", "trial3.rds"),
    dictionary_file = file.path(project_dir, "RealAna", "Trial3", "trial3_dictionary.csv"),
    output_dir = file.path(project_dir, "RealAna", "Trial3", "output"),
    treatment = "Treatment",
    intervention = "Alert",
    control = "Usual Care",
    composite = "YP_composite_event",
    outcomes = c(
      death = "YS_death_event",
      dialysis = "YS_dialysis_event",
      aki = "YS_aki_progression_event"
    ),
    candidates = data.frame(
      Clinical_covariate = c(
        "Hospital / study site", "Age (continuous)", "Age category (fallback)",
        "Sex", "ICU status", "Baseline serum creatinine", "Baseline AKI stage",
        "Elixhauser comorbidity score", "Diabetes", "Congestive heart failure"
      ),
      Variable = c(
        "X_hospital_0m", "X_age_0m", "X_age_over_90_0m",
        "X_sex_0m", "X_icu_0m", "X_creat_0m", NA,
        "X_elx_score_0m", "X_diabetes_0m", "X_chf_0m"
      ),
      Candidate_type = c(
        "core", "age_continuous", "age_category", "core", "core", "core",
        "unavailable", "core", "core", "core"
      ),
      Categorical = c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, NA, FALSE, TRUE, TRUE),
      stringsAsFactors = FALSE
    ),
    baseline_aki_note = "No true baseline AKI-stage variable is present in the dictionary."
  )
)

dictionary_row <- function(dictionary, variable, trial_label) {
  row <- dictionary[dictionary$variable_name == variable, , drop = FALSE]
  assert_true(
    nrow(row) == 1L,
    paste(trial_label, "dictionary mapping is missing or non-unique for", variable)
  )
  row
}

as_binary01 <- function(x, variable, trial_label) {
  values <- as.character(x)
  observed <- sort(unique(values[!is.na(values)]))
  assert_true(
    identical(observed, c("0", "1")),
    paste(trial_label, variable, "must have exactly the nonmissing values 0 and 1.")
  )
  as.integer(values)
}

prepare_trial <- function(config, trial_label) {
  assert_true(file.exists(config$data_file), paste("Missing RDS:", config$data_file))
  assert_true(
    file.exists(config$dictionary_file),
    paste("Missing dictionary:", config$dictionary_file)
  )
  dat <- as.data.frame(readRDS(config$data_file), check.names = FALSE)
  dictionary <- utils::read.csv(config$dictionary_file, check.names = FALSE)
  N <- nrow(dat)
  dat$participant_id <- seq_len(N)

  required <- c(config$treatment, config$composite, unname(config$outcomes))
  missing_required <- setdiff(required, names(dat))
  assert_true(
    length(missing_required) == 0L,
    paste(trial_label, "required variable(s) missing:", paste(missing_required, collapse = ", "))
  )

  treatment_dictionary <- dictionary_row(dictionary, config$treatment, trial_label)
  assert_true(
    treatment_dictionary$variable_role[[1L]] == "Treatment assignment",
    paste(trial_label, "treatment field is not dictionary-labeled as treatment assignment.")
  )
  treatment_values <- unique(as.character(dat[[config$treatment]]))
  assert_true(
    setequal(treatment_values, c(config$control, config$intervention)),
    paste(trial_label, "unexpected treatment levels:", paste(treatment_values, collapse = ", "))
  )
  assert_true(!anyNA(dat[[config$treatment]]), paste(trial_label, "treatment contains NA."))
  dat$A <- as.integer(as.character(dat[[config$treatment]]) == config$intervention)

  for (outcome_name in names(config$outcomes)) {
    variable <- config$outcomes[[outcome_name]]
    row <- dictionary_row(dictionary, variable, trial_label)
    assert_true(
      row$variable_role[[1L]] %in% c("Primary outcome", "Secondary outcome"),
      paste(trial_label, variable, "is not dictionary-labeled as an outcome.")
    )
    dat[[paste0("outcome_", outcome_name)]] <- as_binary01(
      dat[[variable]],
      variable,
      trial_label
    )
    assert_true(
      !anyNA(dat[[paste0("outcome_", outcome_name)]]),
      paste(trial_label, variable, "contains missing outcome values.")
    )
  }
  composite <- as_binary01(dat[[config$composite]], config$composite, trial_label)
  assert_true(
    !anyNA(composite),
    paste(trial_label, config$composite, "contains missing outcome values.")
  )
  component_union <- pmax(
    dat$outcome_death,
    dat$outcome_dialysis,
    dat$outcome_aki
  )
  assert_true(
    identical(composite, component_union),
    paste(
      trial_label,
      "component indicators do not reproduce the trial-defined 14-day composite;",
      "outcome coding is ambiguous."
    )
  )

  candidate_audit <- config$candidates
  candidate_audit$Dictionary_role <- NA_character_
  candidate_audit$Dictionary_description <- NA_character_
  candidate_audit$Missing <- NA_integer_
  candidate_audit$Selected <- FALSE
  candidate_audit$Reason <- ""

  available <- which(!is.na(candidate_audit$Variable))
  for (i in available) {
    variable <- candidate_audit$Variable[[i]]
    assert_true(variable %in% names(dat), paste(trial_label, "candidate missing:", variable))
    row <- dictionary_row(dictionary, variable, trial_label)
    candidate_audit$Dictionary_role[[i]] <- row$variable_role[[1L]]
    candidate_audit$Dictionary_description[[i]] <- row$short_explanation[[1L]]
    candidate_audit$Missing[[i]] <- sum(is.na(dat[[variable]]))
    assert_true(
      row$variable_role[[1L]] == "Baseline covariate",
      paste(trial_label, variable, "is not verified as a baseline covariate.")
    )
  }

  core_rows <- which(candidate_audit$Candidate_type == "core")
  candidate_audit$Selected[core_rows] <- candidate_audit$Missing[core_rows] == 0L
  candidate_audit$Reason[core_rows] <- ifelse(
    candidate_audit$Selected[core_rows],
    "Included: verified baseline and fully observed.",
    paste0(
      "Excluded: ",
      candidate_audit$Missing[core_rows],
      " missing value(s); no imputation or subject deletion."
    )
  )

  age_continuous_row <- which(candidate_audit$Candidate_type == "age_continuous")
  age_category_row <- which(candidate_audit$Candidate_type == "age_category")
  if (candidate_audit$Missing[[age_continuous_row]] == 0L) {
    candidate_audit$Selected[[age_continuous_row]] <- TRUE
    candidate_audit$Reason[[age_continuous_row]] <- "Included: continuous baseline age is complete."
    candidate_audit$Reason[[age_category_row]] <- paste(
      "Not included because complete continuous age was available."
    )
  } else if (candidate_audit$Missing[[age_category_row]] == 0L) {
    candidate_audit$Selected[[age_category_row]] <- TRUE
    candidate_audit$Reason[[age_continuous_row]] <- paste0(
      "Excluded: ", candidate_audit$Missing[[age_continuous_row]],
      " missing value(s); no imputation or subject deletion."
    )
    candidate_audit$Reason[[age_category_row]] <- paste(
      "Included as the fully observed prespecified age-category fallback."
    )
  } else {
    candidate_audit$Reason[[age_continuous_row]] <- "Excluded because age is incomplete."
    candidate_audit$Reason[[age_category_row]] <- "Excluded because age category is incomplete."
  }

  unavailable_row <- which(candidate_audit$Candidate_type == "unavailable")
  candidate_audit$Reason[unavailable_row] <- config$baseline_aki_note

  selected_names <- candidate_audit$Variable[candidate_audit$Selected]
  assert_true(length(selected_names) > 0L, paste(trial_label, "has no LADJ covariates."))
  assert_true(
    all(vapply(dat[selected_names], function(x) !anyNA(x), logical(1))),
    paste(trial_label, "selected LADJ covariate contains NA.")
  )

  categorical_names <- candidate_audit$Variable[
    candidate_audit$Selected & !is.na(candidate_audit$Categorical) &
      candidate_audit$Categorical
  ]

  mapping <- data.frame(
    Role = c(
      "Treatment", "14-day death", "14-day dialysis", "14-day AKI progression",
      candidate_audit$Clinical_covariate[candidate_audit$Selected]
    ),
    Exact_variable = c(
      config$treatment,
      config$outcomes[["death"]],
      config$outcomes[["dialysis"]],
      config$outcomes[["aki"]],
      selected_names
    ),
    Coding_or_use = c(
      paste0("A=1: ", config$intervention, "; A=0: ", config$control),
      "1=death occurrence; 0=no death",
      "1=dialysis occurrence; 0=no dialysis",
      "1=AKI progression occurrence; 0=no progression",
      rep("Baseline LADJ covariate", length(selected_names))
    ),
    stringsAsFactors = FALSE
  )

  cat("\n", strrep("=", 72L), "\n", sep = "")
  cat(trial_label, ": exact mappings verified before estimation\n")
  cat(strrep("=", 72L), "\n", sep = "")
  print(mapping, row.names = FALSE)
  cat("\nCandidate baseline covariate audit:\n")
  print(
    candidate_audit[, c("Clinical_covariate", "Variable", "Missing", "Selected", "Reason")],
    row.names = FALSE
  )
  cat("\n")

  list(
    dat = dat,
    dictionary = dictionary,
    N = N,
    mapping = mapping,
    candidate_audit = candidate_audit,
    selected_names = selected_names,
    categorical_names = categorical_names
  )
}

make_feature_matrix <- function(dat, selected_names, categorical_names) {
  X_data <- dat[selected_names]
  for (variable in categorical_names) {
    observed <- sort(unique(as.character(X_data[[variable]])))
    X_data[[variable]] <- factor(as.character(X_data[[variable]]), levels = observed)
  }
  X <- stats::model.matrix(~ ., data = X_data)
  assert_true("(Intercept)" %in% colnames(X), "Expected model-matrix intercept is absent.")
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  storage.mode(X) <- "double"
  assert_true(!anyNA(X), "LADJ model matrix contains missing values.")
  assert_true(ncol(X) > 0L, "LADJ model matrix has no features.")
  X
}

make_pairwise_comparisons <- function(dat) {
  treated_rows <- which(dat$A == 1L)
  control_rows <- which(dat$A == 0L)
  # Binary weighting is exactly equivalent to lexicographic comparison of
  # (death, dialysis, AKI progression), with lower severity preferred.
  severity <- 4L * dat$outcome_death +
    2L * dat$outcome_dialysis +
    dat$outcome_aki
  severity_t <- severity[treated_rows]
  severity_c <- severity[control_rows]
  win <- outer(severity_t, severity_c, "<")
  loss <- outer(severity_t, severity_c, ">")
  tie_count <- length(win) - sum(win) - sum(loss)

  profile <- data.frame(
    death = dat$outcome_death,
    dialysis = dat$outcome_dialysis,
    aki = dat$outcome_aki
  )
  profile_key <- paste(profile$death, profile$dialysis, profile$aki, sep = "")
  all_keys <- sprintf(
    "%d%d%d",
    rep(0:1, each = 4L),
    rep(rep(0:1, each = 2L), 2L),
    rep(0:1, 4L)
  )
  count_t <- table(factor(profile_key[treated_rows], levels = all_keys))
  count_c <- table(factor(profile_key[control_rows], levels = all_keys))
  profile_grid <- expand.grid(death = 0:1, dialysis = 0:1, aki = 0:1)
  profile_grid <- profile_grid[order(profile_grid$death, profile_grid$dialysis, profile_grid$aki), ]

  tier_counts <- c(death = 0, dialysis = 0, aki = 0, tie = 0)
  for (i in seq_len(nrow(profile_grid))) {
    for (j in seq_len(nrow(profile_grid))) {
      pair_count <- as.numeric(count_t[[i]]) * as.numeric(count_c[[j]])
      if (pair_count == 0) next
      if (profile_grid$death[[i]] != profile_grid$death[[j]]) {
        tier_counts[["death"]] <- tier_counts[["death"]] + pair_count
      } else if (profile_grid$dialysis[[i]] != profile_grid$dialysis[[j]]) {
        tier_counts[["dialysis"]] <- tier_counts[["dialysis"]] + pair_count
      } else if (profile_grid$aki[[i]] != profile_grid$aki[[j]]) {
        tier_counts[["aki"]] <- tier_counts[["aki"]] + pair_count
      } else {
        tier_counts[["tie"]] <- tier_counts[["tie"]] + pair_count
      }
    }
  }

  assert_true(
    sum(tier_counts) == length(win),
    "Decision-tier counts do not partition treatment-control pairs."
  )
  assert_true(tier_counts[["tie"]] == tie_count, "Tie count disagrees with profile count.")
  assert_true(all(win + loss <= 1L), "Pairwise win/loss partition failed.")
  list(
    win = win,
    loss = loss,
    treated_rows = treated_rows,
    control_rows = control_rows,
    severity = severity,
    tier_counts = tier_counts,
    tie_count = tie_count
  )
}

construct_scores <- function(pair_object, N, A, e) {
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
  psi_w[pair_object$treated_rows] <- placement$g_w_T / e
  psi_l[pair_object$treated_rows] <- placement$g_l_T / e
  psi_w[pair_object$control_rows] <- placement$g_w_C / (1 - e)
  psi_l[pair_object$control_rows] <- placement$g_l_C / (1 - e)
  list(
    tau = c(w = tau_w, l = tau_l),
    placement = placement,
    psi = cbind(w = psi_w, l = psi_l)
  )
}

fit_ladj <- function(scores, X, A, e, n1, n0) {
  N <- length(A)
  X_center <- sweep(X, 2L, colMeans(X), "-")
  Z <- A - e
  design <- X_center * Z
  gram <- crossprod(design)
  gram_rank <- qr(gram)$rank
  eigenvalues <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
  assert_true(gram_rank == ncol(X), "Singular LADJ Gram matrix.")
  assert_true(
    min(eigenvalues) > .Machine$double.eps * max(eigenvalues) * ncol(X),
    "Numerically singular LADJ Gram matrix."
  )
  beta <- solve(gram, crossprod(design, scores$psi))
  colnames(beta) <- c("beta_w", "beta_l")
  fitted_b <- X_center %*% beta
  psi_adjusted <- scores$psi - fitted_b * Z
  direct_correction <- colMeans(fitted_b * Z)
  tau_adjusted <- scores$tau - direct_correction
  assert_true(
    all(is.finite(tau_adjusted)) &&
      all(tau_adjusted > 0) &&
      sum(tau_adjusted) < 1,
    "Adjusted elementary win/loss probabilities are outside the admissible region."
  )

  imbalance <- colMeans(X[A == 1L, , drop = FALSE]) -
    colMeans(X[A == 0L, , drop = FALSE])
  imbalance_factor <- n1 * n0 / N^2
  centering_identity <- colMeans(X_center * Z)
  imbalance_correction <- as.vector(imbalance_factor * crossprod(imbalance, beta))

  feature_sd <- apply(X, 2L, stats::sd)
  assert_true(all(feature_sd > 0), "Zero-variance LADJ feature.")
  X_scaled <- sweep(X_center, 2L, feature_sd, "/")
  design_scaled <- X_scaled * Z
  gram_scaled <- crossprod(design_scaled)
  scaled_eigenvalues <- eigen(gram_scaled, symmetric = TRUE, only.values = TRUE)$values
  beta_scaled <- solve(gram_scaled, crossprod(design_scaled, scores$psi))
  fitted_scaled <- X_scaled %*% beta_scaled
  adjusted_covariance <- elementary_covariance(psi_adjusted)

  list(
    tau = tau_adjusted,
    psi = psi_adjusted,
    beta = beta,
    diagnostics = data.frame(
      Number_features = ncol(X),
      Gram_rank = gram_rank,
      Smallest_eigenvalue = min(eigenvalues),
      Largest_eigenvalue = max(eigenvalues),
      Condition_number = max(eigenvalues) / min(eigenvalues),
      Scaled_condition_number = max(scaled_eigenvalues) / min(scaled_eigenvalues),
      Max_centering_identity_error = max(abs(
        centering_identity - imbalance_factor * imbalance
      )),
      Max_tau_formula_error = max(abs(direct_correction - imbalance_correction)),
      Max_scaled_fitted_value_error = max(abs(fitted_b - fitted_scaled)),
      Max_projection_normal_equation_error = max(abs(crossprod(design, psi_adjusted))),
      Min_adjusted_covariance_eigenvalue = min(
        eigen(adjusted_covariance, symmetric = TRUE, only.values = TRUE)$values
      ),
      stringsAsFactors = FALSE
    )
  )
}

elementary_covariance <- function(psi) {
  psi <- as.matrix(psi)
  centered <- sweep(psi, 2L, colMeans(psi), "-")
  crossprod(centered) / nrow(psi)
}

make_inference_rows <- function(trial_label, method, tau, psi_components) {
  N <- nrow(psi_components)
  tau_w <- tau[["w"]]
  tau_l <- tau[["l"]]
  tau_tie <- 1 - tau_w - tau_l
  assert_true(tau_w > 0 && tau_l > 0, "Win/loss probabilities must be positive.")
  delta <- tau_w - tau_l
  assert_true(abs(delta) < 1, "Absolute NB must be less than one.")
  sigma <- elementary_covariance(psi_components)
  z <- stats::qnorm(0.975)

  infer <- function(estimand, estimate, gradient, null_value, scale, transform = NULL) {
    variance <- as.numeric(crossprod(gradient, sigma %*% gradient) / N)
    assert_true(is.finite(variance) && variance >= -1e-14, "Invalid IF variance.")
    se <- sqrt(max(variance, 0))
    ci <- estimate + c(-1, 1) * z * se
    test_statistic <- (estimate - null_value) / se
    p_value <- 2 * stats::pnorm(-abs(test_statistic))
    if (!is.null(transform)) {
      natural_estimate <- transform(estimate)
      natural_ci <- transform(ci)
      natural_se <- natural_estimate * se
      return(data.frame(
        Trial = trial_label,
        Estimand = estimand,
        Method = method,
        Estimate = natural_estimate,
        ASE = natural_se,
        CI_lower = natural_ci[[1L]],
        CI_upper = natural_ci[[2L]],
        P_value = p_value,
        Test_statistic = test_statistic,
        Scale = scale,
        Estimated_variance = natural_se^2,
        tau_w = tau_w,
        tau_l = tau_l,
        tau_tie = tau_tie,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      Trial = trial_label,
      Estimand = estimand,
      Method = method,
      Estimate = estimate,
      ASE = se,
      CI_lower = ci[[1L]],
      CI_upper = ci[[2L]],
      P_value = p_value,
      Test_statistic = test_statistic,
      Scale = scale,
      Estimated_variance = se^2,
      tau_w = tau_w,
      tau_l = tau_l,
      tau_tie = tau_tie,
      stringsAsFactors = FALSE
    )
  }

  grad_nb <- c(1, -1)
  lwr <- log(tau_w) - log(tau_l)
  grad_lwr <- c(1 / tau_w, -1 / tau_l)
  lwo <- log1p(delta) - log1p(-delta)
  grad_lwo <- 2 / (1 - delta^2) * c(1, -1)

  rbind(
    infer("NB", delta, grad_nb, 0, "identity"),
    infer("LWR", lwr, grad_lwr, 0, "log"),
    infer("WR", lwr, grad_lwr, 0, "natural", exp),
    infer("LWO", lwo, grad_lwo, 0, "log"),
    infer("WO", lwo, grad_lwo, 0, "natural", exp)
  )
}

add_efficiency <- function(results) {
  results$ASE_ratio_LADJ_to_UNADJ <- 1
  results$Relative_efficiency_vs_UNADJ <- 1
  results$Percent_variance_reduction <- 0
  for (estimand in unique(results$Estimand)) {
    unadj <- results[results$Estimand == estimand & results$Method == "UNADJ", ]
    ladj <- results[results$Estimand == estimand & results$Method == "LADJ", ]
    assert_true(nrow(unadj) == 1L && nrow(ladj) == 1L, "Incomplete efficiency comparison.")
    ratio <- ladj$ASE / unadj$ASE
    relative_efficiency <- unadj$ASE^2 / ladj$ASE^2
    pvr <- 100 * (1 - ladj$ASE^2 / unadj$ASE^2)
    idx <- results$Estimand == estimand & results$Method == "LADJ"
    results$ASE_ratio_LADJ_to_UNADJ[idx] <- ratio
    results$Relative_efficiency_vs_UNADJ[idx] <- relative_efficiency
    results$Percent_variance_reduction[idx] <- pvr
  }
  results
}

make_arm_summary <- function(dat, config) {
  severity <- ifelse(
    dat$outcome_death == 1L,
    "Death",
    ifelse(
      dat$outcome_dialysis == 1L,
      "Alive + dialysis",
      ifelse(
        dat$outcome_aki == 1L,
        "Alive + no dialysis + AKI progression",
        "None"
      )
    )
  )
  category_levels <- c(
    "Death", "Alive + dialysis", "Alive + no dialysis + AKI progression", "None"
  )
  do.call(rbind, lapply(c(config$control, config$intervention), function(arm) {
    keep <- as.character(dat[[config$treatment]]) == arm
    n <- sum(keep)
    category_counts <- table(factor(severity[keep], levels = category_levels))
    data.frame(
      Arm = arm,
      Role = ifelse(arm == config$intervention, "Intervention", "Control"),
      N = n,
      Death_n = sum(dat$outcome_death[keep]),
      Death_percent = 100 * mean(dat$outcome_death[keep]),
      Dialysis_n = sum(dat$outcome_dialysis[keep]),
      Dialysis_percent = 100 * mean(dat$outcome_dialysis[keep]),
      AKI_progression_n = sum(dat$outcome_aki[keep]),
      AKI_progression_percent = 100 * mean(dat$outcome_aki[keep]),
      Severity_death_n = unname(category_counts[["Death"]]),
      Severity_death_percent = 100 * unname(category_counts[["Death"]]) / n,
      Severity_alive_dialysis_n = unname(category_counts[["Alive + dialysis"]]),
      Severity_alive_dialysis_percent =
        100 * unname(category_counts[["Alive + dialysis"]]) / n,
      Severity_alive_no_dialysis_AKI_n =
        unname(category_counts[["Alive + no dialysis + AKI progression"]]),
      Severity_alive_no_dialysis_AKI_percent =
        100 * unname(category_counts[["Alive + no dialysis + AKI progression"]]) / n,
      Severity_none_n = unname(category_counts[["None"]]),
      Severity_none_percent = 100 * unname(category_counts[["None"]]) / n,
      stringsAsFactors = FALSE
    )
  }))
}

make_balance <- function(X, A) {
  do.call(rbind, lapply(seq_len(ncol(X)), function(j) {
    x <- X[, j]
    n1 <- sum(A == 1L)
    n0 <- sum(A == 0L)
    pooled_sd <- sqrt(
      ((n1 - 1) * stats::var(x[A == 1L]) + (n0 - 1) * stats::var(x[A == 0L])) /
        (n1 + n0 - 2)
    )
    data.frame(
      Feature = colnames(X)[[j]],
      Intervention_mean = mean(x[A == 1L]),
      Control_mean = mean(x[A == 0L]),
      Raw_mean_difference = mean(x[A == 1L]) - mean(x[A == 0L]),
      Pooled_within_arm_SD = pooled_sd,
      Standardized_mean_difference = if (pooled_sd > 0) {
        (mean(x[A == 1L]) - mean(x[A == 0L])) / pooled_sd
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
}

analyze_trial <- function(config, trial_label, prepared) {
  dat <- prepared$dat
  N <- prepared$N
  n1 <- sum(dat$A == 1L)
  n0 <- sum(dat$A == 0L)
  assert_true(n1 > 0L && n0 > 0L, paste(trial_label, "does not contain both arms."))

  X <- make_feature_matrix(
    dat,
    prepared$selected_names,
    prepared$categorical_names
  )
  pair_object <- make_pairwise_comparisons(dat)
  scores <- construct_scores(pair_object, N, dat$A, known_e)
  ladj <- fit_ladj(scores, X, dat$A, known_e, n1, n0)

  unadj_results <- make_inference_rows(
    trial_label,
    "UNADJ",
    scores$tau,
    scores$psi
  )
  ladj_results <- make_inference_rows(
    trial_label,
    "LADJ",
    ladj$tau,
    ladj$psi
  )
  results <- add_efficiency(rbind(unadj_results, ladj_results))
  estimand_order <- c("NB", "LWR", "WR", "LWO", "WO")
  results <- results[
    order(match(results$Estimand, estimand_order), match(results$Method, c("UNADJ", "LADJ"))),
  ]

  N_pairs <- n1 * n0
  pairwise_summary <- data.frame(
    Trial = trial_label,
    Horizon_days = analysis_horizon_days,
    N = N,
    n1 = n1,
    n0 = n0,
    N_pairs = N_pairs,
    Treatment_wins = sum(pair_object$win),
    Treatment_losses = sum(pair_object$loss),
    Final_ties = pair_object$tie_count,
    tau_w = scores$tau[["w"]],
    tau_l = scores$tau[["l"]],
    tau_tie = 1 - sum(scores$tau),
    Decided_at_death_n = pair_object$tier_counts[["death"]],
    Decided_at_death_percent = 100 * pair_object$tier_counts[["death"]] / N_pairs,
    Decided_at_dialysis_n = pair_object$tier_counts[["dialysis"]],
    Decided_at_dialysis_percent = 100 * pair_object$tier_counts[["dialysis"]] / N_pairs,
    Decided_at_AKI_progression_n = pair_object$tier_counts[["aki"]],
    Decided_at_AKI_progression_percent = 100 * pair_object$tier_counts[["aki"]] / N_pairs,
    Final_ties_percent = 100 * pair_object$tier_counts[["tie"]] / N_pairs,
    stringsAsFactors = FALSE
  )
  assert_true(
    pairwise_summary$Treatment_wins +
      pairwise_summary$Treatment_losses +
      pairwise_summary$Final_ties == N_pairs,
    paste(trial_label, "pair counts do not sum to n1*n0.")
  )

  arm_summary <- make_arm_summary(dat, config)
  covariate_audit <- prepared$candidate_audit
  mapping <- prepared$mapping
  outcome_missingness <- data.frame(
    Outcome = c("14-day death", "14-day dialysis", "14-day AKI progression"),
    Variable = unname(config$outcomes),
    Missing = vapply(unname(config$outcomes), function(v) sum(is.na(dat[[v]])), integer(1)),
    stringsAsFactors = FALSE
  )
  population <- data.frame(
    Trial = trial_label,
    Population = c("All randomized", config$control, config$intervention),
    Role = c("Total", "Control", "Intervention"),
    N = c(N, n0, n1),
    Treatment_code_A = c(NA, 0, 1),
    Known_randomization_probability = c(known_e, 1 - known_e, known_e),
    stringsAsFactors = FALSE
  )
  elementary <- data.frame(
    Trial = trial_label,
    Method = c("UNADJ", "LADJ"),
    tau_w = c(scores$tau[["w"]], ladj$tau[["w"]]),
    tau_l = c(scores$tau[["l"]], ladj$tau[["l"]]),
    tau_tie = c(1 - sum(scores$tau), 1 - sum(ladj$tau)),
    stringsAsFactors = FALSE
  )
  coefficients <- data.frame(
    Feature = rownames(ladj$beta),
    beta_w = ladj$beta[, "beta_w"],
    beta_l = ladj$beta[, "beta_l"],
    stringsAsFactors = FALSE
  )
  diagnostics <- cbind(Trial = trial_label, ladj$diagnostics)
  balance <- make_balance(X, dat$A)

  feature_to_raw <- vapply(colnames(X), function(feature) {
    hits <- prepared$selected_names[startsWith(feature, prepared$selected_names)]
    assert_true(
      length(hits) > 0L,
      paste("Could not map LADJ design feature to raw covariate:", feature)
    )
    hits[[which.max(nchar(hits))]]
  }, character(1))
  imbalance_factor <- n1 * n0 / N^2
  correction_decomposition <- merge(
    balance,
    coefficients,
    by.x = "Feature",
    by.y = "Feature",
    sort = FALSE
  )
  correction_decomposition$Raw_covariate <- feature_to_raw[
    match(correction_decomposition$Feature, colnames(X))
  ]
  correction_decomposition$Delta_tau_w <- -imbalance_factor *
    correction_decomposition$Raw_mean_difference *
    correction_decomposition$beta_w
  correction_decomposition$Delta_tau_l <- -imbalance_factor *
    correction_decomposition$Raw_mean_difference *
    correction_decomposition$beta_l
  correction_decomposition$First_order_delta_LWR <- (
    correction_decomposition$Delta_tau_w / scores$tau[["w"]] -
      correction_decomposition$Delta_tau_l / scores$tau[["l"]]
  )
  correction_decomposition <- correction_decomposition[
    order(-abs(correction_decomposition$First_order_delta_LWR)),
  ]

  unadj_lwr <- results[
    results$Estimand == "LWR" & results$Method == "UNADJ",
  ]
  leave_one_out <- lapply(c("Full primary", prepared$selected_names), function(label) {
    if (label == "Full primary") {
      fit <- ladj
      omitted <- ""
    } else {
      keep <- feature_to_raw != label
      fit <- fit_ladj(
        scores,
        X[, keep, drop = FALSE],
        dat$A,
        known_e,
        n1,
        n0
      )
      omitted <- label
    }
    lwr <- make_inference_rows(trial_label, "LADJ", fit$tau, fit$psi)
    lwr <- lwr[lwr$Estimand == "LWR", ]
    data.frame(
      Model = if (label == "Full primary") "Full primary LADJ" else paste("Omit", label),
      Omitted_raw_covariate = omitted,
      LWR = lwr$Estimate,
      ASE = lwr$ASE,
      CI_lower = lwr$CI_lower,
      CI_upper = lwr$CI_upper,
      P_value = lwr$P_value,
      Percent_variance_reduction = 100 * (1 - (lwr$ASE / unadj_lwr$ASE)^2),
      tau_w = fit$tau[["w"]],
      tau_l = fit$tau[["l"]],
      tau_tie = 1 - sum(fit$tau),
      stringsAsFactors = FALSE
    )
  })
  leave_one_out <- do.call(rbind, leave_one_out)

  realized_e <- mean(dat$A)
  realized_scores <- construct_scores(pair_object, N, dat$A, realized_e)
  realized_ladj <- fit_ladj(realized_scores, X, dat$A, realized_e, n1, n0)
  realized_lwr <- make_inference_rows(
    trial_label,
    "LADJ",
    realized_ladj$tau,
    realized_ladj$psi
  )
  realized_lwr <- realized_lwr[realized_lwr$Estimand == "LWR", ]
  e_sensitivity <- rbind(
    data.frame(
      Probability_definition = "Known randomization probability",
      e = known_e,
      LWR = results$Estimate[
        results$Estimand == "LWR" & results$Method == "LADJ"
      ],
      ASE = results$ASE[results$Estimand == "LWR" & results$Method == "LADJ"],
      stringsAsFactors = FALSE
    ),
    data.frame(
      Probability_definition = "Realized treatment proportion (diagnostic only)",
      e = realized_e,
      LWR = realized_lwr$Estimate,
      ASE = realized_lwr$ASE,
      stringsAsFactors = FALSE
    )
  )

  max_tau_transform_error <- max(vapply(c("UNADJ", "LADJ"), function(method) {
    method_rows <- results[results$Method == method, ]
    tau_w <- method_rows$tau_w[[1L]]
    tau_l <- method_rows$tau_l[[1L]]
    delta <- tau_w - tau_l
    expected <- c(
      NB = delta,
      LWR = log(tau_w / tau_l),
      WR = tau_w / tau_l,
      LWO = log((1 + delta) / (1 - delta)),
      WO = (1 + delta) / (1 - delta)
    )
    max(abs(method_rows$Estimate[match(names(expected), method_rows$Estimand)] - expected))
  }, numeric(1)))

  qc <- data.frame(
    Check = c(
      "All randomized subjects retained",
      "Treatment has no missing values and both arms are present",
      "All three 14-day outcomes are complete",
      "Outcome indicators are binary 0/1",
      "Component union reproduces trial-defined 14-day composite",
      "Every selected baseline covariate is complete",
      "No post-randomization variable is used for LADJ",
      "Pairwise win/loss/tie partition",
      "LADJ Gram matrix full rank",
      "LADJ projection normal equations",
      "Adjusted elementary covariance is positive semidefinite",
      "Exact transformations from elementary tau"
    ),
    Passed = c(
      nrow(dat) == N,
      !anyNA(dat[[config$treatment]]) && n1 > 0L && n0 > 0L,
      all(outcome_missingness$Missing == 0L),
      all(dat$outcome_death %in% 0:1) &&
        all(dat$outcome_dialysis %in% 0:1) &&
        all(dat$outcome_aki %in% 0:1),
      identical(
        as_binary01(dat[[config$composite]], config$composite, trial_label),
        pmax(dat$outcome_death, dat$outcome_dialysis, dat$outcome_aki)
      ),
      all(vapply(dat[prepared$selected_names], function(x) !anyNA(x), logical(1))),
      all(covariate_audit$Dictionary_role[covariate_audit$Selected] == "Baseline covariate"),
      pairwise_summary$Treatment_wins +
        pairwise_summary$Treatment_losses +
        pairwise_summary$Final_ties == N_pairs,
      diagnostics$Gram_rank == diagnostics$Number_features,
      diagnostics$Max_projection_normal_equation_error < 1e-9,
      diagnostics$Min_adjusted_covariance_eigenvalue > -1e-12,
      max_tau_transform_error < 1e-12
    ),
    Details = c(
      paste("N =", N),
      paste("n1 =", n1, "; n0 =", n0),
      paste(outcome_missingness$Variable, outcome_missingness$Missing, collapse = "; "),
      "Dictionary coding and observed values verified",
      config$composite,
      paste(prepared$selected_names, collapse = ", "),
      "All selected fields have dictionary role Baseline covariate",
      paste("n1*n0 =", N_pairs),
      paste("rank =", diagnostics$Gram_rank, "of", diagnostics$Number_features),
      paste(
        "maximum normal-equation residual =",
        format(diagnostics$Max_projection_normal_equation_error, scientific = TRUE)
      ),
      paste(
        "minimum covariance eigenvalue =",
        format(diagnostics$Min_adjusted_covariance_eigenvalue, scientific = TRUE)
      ),
      paste("maximum error =", format(max_tau_transform_error, scientific = TRUE))
    ),
    stringsAsFactors = FALSE
  )
  assert_true(all(qc$Passed), paste(trial_label, "quality-control check failed."))

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- paste0("trial", config$trial_id, "_")
  write_csv(mapping, file.path(config$output_dir, paste0(prefix, "variable_mapping.csv")))
  write_csv(population, file.path(config$output_dir, paste0(prefix, "population.csv")))
  write_csv(
    outcome_missingness,
    file.path(config$output_dir, paste0(prefix, "outcome_missingness.csv"))
  )
  write_csv(
    arm_summary,
    file.path(config$output_dir, paste0(prefix, "tableA_arm_outcomes.csv"))
  )
  write_csv(
    covariate_audit,
    file.path(config$output_dir, paste0(prefix, "tableB_covariates.csv"))
  )
  write_csv(results, file.path(config$output_dir, paste0(prefix, "tableC_results.csv")))
  write_csv(
    pairwise_summary,
    file.path(config$output_dir, paste0(prefix, "pairwise_summary.csv"))
  )
  write_csv(
    elementary,
    file.path(config$output_dir, paste0(prefix, "elementary_estimates.csv"))
  )
  write_csv(
    coefficients,
    file.path(config$output_dir, paste0(prefix, "ladj_coefficients.csv"))
  )
  write_csv(
    diagnostics,
    file.path(config$output_dir, paste0(prefix, "ladj_diagnostics.csv"))
  )
  write_csv(balance, file.path(config$output_dir, paste0(prefix, "baseline_balance.csv")))
  write_csv(
    correction_decomposition,
    file.path(config$output_dir, paste0(prefix, "ladj_correction_decomposition.csv"))
  )
  write_csv(
    leave_one_out,
    file.path(config$output_dir, paste0(prefix, "ladj_leave_one_covariate_out.csv"))
  )
  write_csv(
    e_sensitivity,
    file.path(config$output_dir, paste0(prefix, "allocation_probability_diagnostic.csv"))
  )
  write_csv(qc, file.path(config$output_dir, paste0(prefix, "quality_control.csv")))

  report_results <- results[, c(
    "Trial", "Estimand", "Method", "Estimate", "ASE", "CI_lower", "CI_upper",
    "P_value", "ASE_ratio_LADJ_to_UNADJ", "Relative_efficiency_vs_UNADJ",
    "Percent_variance_reduction"
  )]
  report_lines <- c(
    paste0("# ", trial_label, ": 14-day hierarchical win analysis"),
    "",
    "## Verified mappings",
    "",
    markdown_table(mapping),
    "",
    "## Analysis population",
    "",
    markdown_table(population),
    "",
    "## Outcome missingness",
    "",
    markdown_table(outcome_missingness),
    "",
    "## Table A: outcomes by randomized arm",
    "",
    markdown_table(arm_summary),
    "",
    "## Unadjusted pairwise comparison",
    "",
    markdown_table(pairwise_summary),
    "",
    "## Table B: prespecified baseline covariate audit",
    "",
    markdown_table(covariate_audit[, c(
      "Clinical_covariate", "Variable", "Missing", "Selected", "Reason"
    )]),
    "",
    paste0(
      "**Final LADJ raw covariates:** ",
      paste(prepared$selected_names, collapse = ", "),
      "."
    ),
    "",
    paste0(
      "**Final LADJ design columns after reference dummy coding:** ",
      paste(colnames(X), collapse = ", "),
      "."
    ),
    "",
    "## Table C: unadjusted and LADJ estimates",
    "",
    markdown_table(report_results, digits = 6L),
    "",
    "## Elementary win/loss estimates",
    "",
    markdown_table(elementary, digits = 6L),
    "",
    "## LADJ diagnostics",
    "",
    markdown_table(diagnostics, digits = 6L),
    "",
    "## LADJ point-estimate correction diagnostic",
    "",
    paste(
      "The feature-level LWR contributions below are a first-order decomposition of the",
      "exact elementary-probability correction. They are diagnostic and do not redefine",
      "the prespecified primary model."
    ),
    "",
    markdown_table(correction_decomposition[, c(
      "Feature", "Raw_covariate", "Standardized_mean_difference",
      "Delta_tau_w", "Delta_tau_l", "First_order_delta_LWR"
    )], digits = 6L),
    "",
    "## Leave-one-covariate-out diagnostic",
    "",
    paste(
      "These post-hoc fits diagnose the source of the adjustment and are not",
      "outcome-driven alternatives to the primary LADJ analysis."
    ),
    "",
    markdown_table(leave_one_out, digits = 6L),
    "",
    "## Allocation-probability diagnostic",
    "",
    markdown_table(e_sensitivity, digits = 6L),
    "",
    "## Quality control",
    "",
    markdown_table(qc)
  )
  writeLines(
    report_lines,
    file.path(config$output_dir, paste0(prefix, "analysis_summary.md"))
  )

  cat(trial_label, "analysis complete. Output:", config$output_dir, "\n")
  list(
    trial = trial_label,
    config = config,
    results = results,
    pairwise = pairwise_summary,
    elementary = elementary,
    covariates = covariate_audit,
    diagnostics = diagnostics,
    arm_summary = arm_summary
  )
}

trial_label <- "AKI alert trial"
prepared <- prepare_trial(trial_configs$Trial3, trial_label)
analysis <- analyze_trial(trial_configs$Trial3, trial_label, prepared)

cat("\nAKI alert trial analysis completed.\n")
