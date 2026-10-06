required_sim_packages <- function() {
  c("data.table", "future.apply", "future", "splines")
}

check_sim_packages <- function() {
  missing <- required_sim_packages()[!vapply(required_sim_packages(), requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing required R package(s): ",
      paste(missing, collapse = ", "),
      ". Install them before running the simulation.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

default_dgm_params <- function() {
  S <- 48
  weibull_k <- 1.2
  alternative_effect_scale <- 0.70

  list(
    S = S,
    # Weibull hazard is lambda0 * t^(k - 1) * exp(zeta_D1).
    weibull_k = weibull_k,
    # Calibrated at zeta_D1 = 0 to give 20% death by S:
    # exp{-lambda0 * S^k / k} = 0.80.
    lambda0 = weibull_k * (-log(0.80)) / (S^weibull_k),
    theta_D10 = 0,
    theta_D20 = 36,
    theta_D30 = log(60),
    death_epsilon = 0,
    # Treat ALSFRS-S differences within 2 points as clinically tied.
    alsfrs_epsilon = 2,
    # NfL is compared on -log(NfL); log(1.25) means a within-25% ratio is tied.
    nfl_epsilon = log(1.25),
    sigma_D2 = 5,
    sigma_D3 = 0.45,
    alpha_null = c(D1 = 0, D2 = 0, D3 = 0),
    # Moderate alternative used across all sample-size configurations.
    alternative_effect_scale = alternative_effect_scale,
    alpha_alt = alternative_effect_scale * c(D1 = log(0.875), D2 = 1.15, D3 = log(0.935)),
    prognosis_scale = c(weak = 0.35, moderate = 0.70, strong = 1.10),
    heterogeneity_scale = c(homogeneous = 0, heterogeneous = 1),
    beta_D1_linear_base = c(X1 = 0.45, X2 = -0.35, X3 = 0.55),
    beta_D2_linear_base = c(X1 = -3.0, X2 = 2.0, X3 = -4.0),
    beta_D3_linear_base = c(X1 = 0.30, X2 = -0.22, X3 = 0.35),
    beta_D1_endpoint_specific_base = c(X1 = 0.45, X2 = -0.35, X3 = 0.00),
    beta_D2_endpoint_specific_base = c(X1 = 0.00, X2 = 2.00, X3 = -4.00),
    beta_D3_endpoint_specific_base = c(X1 = 0.00, X2 = 0.00, X3 = 0.35),
    beta_D1_nonlinear_base = c(X1_sq = 0.25, X1_X3 = 0.35),
    beta_D2_nonlinear_base = c(X1_sq = -1.8, X1_X3 = -2.2),
    beta_D3_nonlinear_base = c(X1_sq = 0.18, X1_X3 = 0.22),
    gamma_D1_base = alternative_effect_scale * c(X1 = -0.08, X2 = 0.05, X3 = 0.06),
    gamma_D2_base = alternative_effect_scale * c(X1 = 0.25, X2 = -0.15, X3 = 0.20),
    gamma_D3_base = alternative_effect_scale * c(X1 = -0.035, X2 = 0.025, X3 = -0.025)
  )
}

normalize_scenario_grid <- function(grid) {
  grid <- data.table::as.data.table(grid)

  if (!"outcome_shape" %in% names(grid)) {
    if ("specification" %in% names(grid)) {
      grid[, outcome_shape := ifelse(
        specification %in% c("misspecified", "nonlinear"),
        "nonlinear",
        "linear"
      )]
    } else {
      grid[, outcome_shape := "linear"]
    }
  }

  if (!"effect_heterogeneity" %in% names(grid)) {
    grid[, effect_heterogeneity := "homogeneous"]
  }
  if (!"endpoint_prognosis_pattern" %in% names(grid)) {
    grid[, endpoint_prognosis_pattern := "shared"]
  }
  if (!"n_noise_covariates" %in% names(grid)) {
    grid[, n_noise_covariates := 0L]
  }
  if (!"adjustment_covariate_pattern" %in% names(grid)) {
    grid[, adjustment_covariate_pattern := "all_available"]
  }
  grid[, n_noise_covariates := as.integer(n_noise_covariates)]

  bad_shape <- setdiff(unique(grid$outcome_shape), c("linear", "nonlinear"))
  if (length(bad_shape) > 0L) {
    stop("Unknown outcome_shape value(s): ", paste(bad_shape, collapse = ", "), call. = FALSE)
  }

  bad_heterogeneity <- setdiff(unique(grid$effect_heterogeneity), c("homogeneous", "heterogeneous"))
  if (length(bad_heterogeneity) > 0L) {
    stop(
      "Unknown effect_heterogeneity value(s): ",
      paste(bad_heterogeneity, collapse = ", "),
      call. = FALSE
    )
  }

  bad_endpoint_pattern <- setdiff(unique(grid$endpoint_prognosis_pattern), c("shared", "endpoint_specific"))
  if (length(bad_endpoint_pattern) > 0L) {
    stop(
      "Unknown endpoint_prognosis_pattern value(s): ",
      paste(bad_endpoint_pattern, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(grid$n_noise_covariates) || any(!grid$n_noise_covariates %in% 0:2)) {
    stop("n_noise_covariates must be 0, 1, or 2.", call. = FALSE)
  }

  bad_adjustment_pattern <- setdiff(
    unique(grid$adjustment_covariate_pattern),
    c("all_available", "omit_X3_add_X4", "omit_X3")
  )
  if (length(bad_adjustment_pattern) > 0L) {
    stop(
      "Unknown adjustment_covariate_pattern value(s): ",
      paste(bad_adjustment_pattern, collapse = ", "),
      call. = FALSE
    )
  }
  if (any(
    grid$adjustment_covariate_pattern == "omit_X3_add_X4" &
      grid$n_noise_covariates < 1L
  )) {
    stop("omit_X3_add_X4 requires at least one noise covariate.", call. = FALSE)
  }

  if (!"scenario_id" %in% names(grid)) {
    grid[, scenario_id := .I]
  }
  if (!"rng_scenario_id" %in% names(grid)) {
    grid[, rng_scenario_id := scenario_id]
  } else {
    grid[is.na(rng_scenario_id), rng_scenario_id := scenario_id]
  }
  grid[, rng_scenario_id := as.integer(rng_scenario_id)]
  if (anyNA(grid$rng_scenario_id) || any(grid$rng_scenario_id < 1L)) {
    stop("rng_scenario_id must contain positive integers.", call. = FALSE)
  }

  grid[]
}

make_full_factorial_scenario_grid <- function() {
  grid <- expand.grid(
    N = c(200L, 600L, 1000L),
    e = c(0.5, 2 / 3),
    effect = c("null", "alternative"),
    prognosis = c("weak", "moderate", "strong"),
    outcome_shape = c("linear", "nonlinear"),
    effect_heterogeneity = c("homogeneous", "heterogeneous"),
    endpoint_prognosis_pattern = c("shared", "endpoint_specific"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  grid <- data.table::as.data.table(grid)
  data.table::setorder(grid, effect, prognosis, outcome_shape, effect_heterogeneity, endpoint_prognosis_pattern, e, N)
  grid[, scenario_id := .I]
  grid[]
}

make_one_factor_scenario_grid <- function(run_full = FALSE) {
  baseline <- data.table::data.table(
    N = 1000L,
    e = 0.5,
    prognosis = "moderate",
    outcome_shape = "linear",
    effect_heterogeneity = "homogeneous",
    endpoint_prognosis_pattern = "shared",
    n_noise_covariates = 0L,
    adjustment_covariate_pattern = "all_available",
    scenario_change = "baseline"
  )

  scenarios <- data.table::rbindlist(list(
    baseline,
    data.table::copy(baseline)[, `:=`(N = 500L, scenario_change = "N_500")],
    data.table::copy(baseline)[, `:=`(e = 2 / 3, scenario_change = "allocation_1_to_2")],
    data.table::copy(baseline)[, `:=`(prognosis = "weak", scenario_change = "weak_prognosis")],
    data.table::copy(baseline)[, `:=`(prognosis = "strong", scenario_change = "strong_prognosis")],
    data.table::copy(baseline)[, `:=`(
      effect_heterogeneity = "heterogeneous",
      scenario_change = "heterogeneous_effect"
    )],
    data.table::copy(baseline)[, `:=`(
      endpoint_prognosis_pattern = "endpoint_specific",
      scenario_change = "endpoint_specific_prognosis"
    )],
    data.table::copy(baseline)[, `:=`(
      n_noise_covariates = 2L,
      scenario_change = "two_noise_covariates"
    )],
    data.table::copy(baseline)[, `:=`(
      adjustment_covariate_pattern = "omit_X3",
      scenario_change = "omit_X3"
    )]
  ), use.names = TRUE)

  grid <- data.table::rbindlist(lapply(c("null", "alternative"), function(effect_value) {
    out <- data.table::copy(scenarios)
    out[, effect := effect_value]
    out
  }), use.names = TRUE)
  grid[, scenario_id := .I]
  grid[, rng_scenario_id := scenario_id]
  grid[scenario_change %in% c("two_noise_covariates", "omit_X3"),
       rng_scenario_id := ifelse(effect == "null", 1L, 10L)]
  grid[]
}

make_scenario_grid <- function(run_full = FALSE,
                               scenario_design = c("one_factor", "full_factorial")) {
  scenario_design <- match.arg(scenario_design)
  if (identical(scenario_design, "full_factorial")) {
    make_full_factorial_scenario_grid()
  } else {
    make_one_factor_scenario_grid(run_full)
  }
}

simulation_settings <- function(run_full = FALSE) {
  if (run_full) {
    list(
      n_rep = 20000L,
      truth_pairs = 10000000L,
      truth_chunk_size = 100000L,
      sbw_nonnegative = FALSE,
      bootstrap_B = 0L
    )
  } else {
    list(
      n_rep = 20L,
      truth_pairs = 20000L,
      truth_chunk_size = 5000L,
      sbw_nonnegative = FALSE,
      bootstrap_B = 0L
    )
  }
}

adjustment_covariate_names <- function(dat) {
  covariate_names <- attr(dat, "adjustment_covariates", exact = TRUE)
  if (is.null(covariate_names)) {
    covariate_names <- intersect(paste0("X", seq_len(5L)), names(dat))
  }
  missing <- setdiff(covariate_names, names(dat))
  if (length(missing) > 0L) {
    stop("Adjustment covariate(s) missing from data: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  covariate_names
}

adjustment_features <- function(dat) {
  covariate_names <- adjustment_covariate_names(dat)
  X <- as.matrix(dat[, ..covariate_names])
  sweep(X, 2L, colMeans(X), FUN = "-")
}
