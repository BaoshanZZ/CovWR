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
    # A moderate alternative keeps N = 1000 power below saturation in the main grid.
    alpha_alt = c(D1 = log(0.875), D2 = 1.15, D3 = log(0.935)),
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
    gamma_D1_base = c(X1 = -0.08, X2 = 0.05, X3 = 0.06),
    gamma_D2_base = c(X1 = 0.25, X2 = -0.15, X3 = 0.20),
    gamma_D3_base = c(X1 = -0.035, X2 = 0.025, X3 = -0.025),
    min_tau = 1e-8
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

  if (!"scenario_id" %in% names(grid)) {
    grid[, scenario_id := .I]
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
    endpoint_prognosis_pattern = "shared"
  )

  if (run_full) {
    scenarios <- data.table::rbindlist(list(
      baseline,
      baseline[, .(N = c(500L, 1000L), e, prognosis, outcome_shape, effect_heterogeneity, endpoint_prognosis_pattern)],
      baseline[, .(N, e, prognosis = c("weak", "moderate", "strong"), outcome_shape, effect_heterogeneity, endpoint_prognosis_pattern)],
      baseline[, .(N, e, prognosis, outcome_shape = c("linear", "nonlinear"), effect_heterogeneity, endpoint_prognosis_pattern)],
      baseline[, .(N, e, prognosis, outcome_shape, effect_heterogeneity = c("homogeneous", "heterogeneous"), endpoint_prognosis_pattern)],
      baseline[, .(N, e, prognosis, outcome_shape, effect_heterogeneity, endpoint_prognosis_pattern = c("shared", "endpoint_specific"))]
    ), use.names = TRUE)
    effects <- c("null", "alternative")
  } else {
    scenarios <- data.table::rbindlist(list(
      baseline,
      baseline[, .(N = 500L, e, prognosis, outcome_shape, effect_heterogeneity, endpoint_prognosis_pattern)],
      baseline[, .(N, e, prognosis = c("weak", "strong"), outcome_shape, effect_heterogeneity, endpoint_prognosis_pattern)],
      baseline[, .(N, e, prognosis, outcome_shape = "nonlinear", effect_heterogeneity, endpoint_prognosis_pattern)],
      baseline[, .(N, e, prognosis, outcome_shape, effect_heterogeneity = "heterogeneous", endpoint_prognosis_pattern)],
      baseline[, .(N, e, prognosis, outcome_shape, effect_heterogeneity, endpoint_prognosis_pattern = "endpoint_specific")]
    ), use.names = TRUE)
    effects <- c("null", "alternative")
  }

  scenarios <- unique(scenarios)
  grid <- scenarios[, .(effect = effects), by = .(
    N, e, prognosis, outcome_shape, effect_heterogeneity, endpoint_prognosis_pattern
  )]
  data.table::setorder(grid, effect, prognosis, outcome_shape, effect_heterogeneity, endpoint_prognosis_pattern, e, N)
  grid[, scenario_id := .I]
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
      n_rep = 5000L,
      truth_pairs = 1000000L,
      truth_chunk_size = 100000L,
      sbw_nonnegative = FALSE
    )
  } else {
    list(
      n_rep = 20L,
      truth_pairs = 20000L,
      truth_chunk_size = 5000L,
      sbw_nonnegative = FALSE
    )
  }
}

adjustment_features <- function(dat) {
  cbind(
    intercept = 1,
    X1 = dat$X1,
    X2 = dat$X2,
    X3 = dat$X3
  )
}
