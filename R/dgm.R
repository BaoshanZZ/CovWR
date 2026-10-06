get_alpha_effect <- function(effect, params = default_dgm_params()) {
  effect <- match.arg(effect, c("null", "alternative"))
  if (identical(effect, "null")) {
    params$alpha_null
  } else {
    params$alpha_alt
  }
}

get_prognosis_scale <- function(prognosis, params = default_dgm_params()) {
  prognosis <- match.arg(prognosis, names(params$prognosis_scale))
  unname(params$prognosis_scale[[prognosis]])
}

get_linear_endpoint_betas <- function(endpoint_prognosis_pattern,
                                      params = default_dgm_params()) {
  endpoint_prognosis_pattern <- match.arg(endpoint_prognosis_pattern, c("shared", "endpoint_specific"))
  if (identical(endpoint_prognosis_pattern, "endpoint_specific")) {
    list(
      D1 = params$beta_D1_endpoint_specific_base,
      D2 = params$beta_D2_endpoint_specific_base,
      D3 = params$beta_D3_endpoint_specific_base
    )
  } else {
    list(
      D1 = params$beta_D1_linear_base,
      D2 = params$beta_D2_linear_base,
      D3 = params$beta_D3_linear_base
    )
  }
}

make_nonlinear_terms <- function(X1, X3) {
  cbind(
    X1_sq = X1^2 - 1,
    X1_X3 = X1 * X3
  )
}

get_heterogeneity_scale <- function(effect,
                                    effect_heterogeneity,
                                    params = default_dgm_params()) {
  effect <- match.arg(effect, c("null", "alternative"))
  effect_heterogeneity <- match.arg(effect_heterogeneity, names(params$heterogeneity_scale))
  if (identical(effect, "null")) {
    return(0)
  }
  unname(params$heterogeneity_scale[[effect_heterogeneity]])
}

map_legacy_specification <- function(specification) {
  specification <- match.arg(specification, c("correct", "misspecified", "linear", "nonlinear"))
  if (specification %in% c("misspecified", "nonlinear")) {
    "nonlinear"
  } else {
    "linear"
  }
}

generate_subject_data <- function(N,
                                  A = NULL,
                                  e = 0.5,
                                  effect = "null",
                                  prognosis = "moderate",
                                  outcome_shape = "linear",
                                  effect_heterogeneity = "homogeneous",
                                  endpoint_prognosis_pattern = "shared",
                                  n_noise_covariates = 0L,
                                  adjustment_covariate_pattern = "all_available",
                                  specification = NULL,
                                  params = default_dgm_params()) {
  effect <- match.arg(effect, c("null", "alternative"))
  if (!is.null(specification)) {
    outcome_shape <- map_legacy_specification(specification)
  }
  outcome_shape <- match.arg(outcome_shape, c("linear", "nonlinear"))
  effect_heterogeneity <- match.arg(effect_heterogeneity, c("homogeneous", "heterogeneous"))
  endpoint_prognosis_pattern <- match.arg(endpoint_prognosis_pattern, c("shared", "endpoint_specific"))
  n_noise_covariates <- as.integer(n_noise_covariates)
  if (is.na(n_noise_covariates) || !n_noise_covariates %in% 0:2) {
    stop("n_noise_covariates must be 0, 1, or 2.", call. = FALSE)
  }
  adjustment_covariate_pattern <- match.arg(
    adjustment_covariate_pattern,
    c("all_available", "omit_X3_add_X4", "omit_X3")
  )
  if (identical(adjustment_covariate_pattern, "omit_X3_add_X4") && n_noise_covariates < 1L) {
    stop("omit_X3_add_X4 requires at least one noise covariate.", call. = FALSE)
  }

  if (is.null(A)) {
    A <- stats::rbinom(N, size = 1L, prob = e)
  } else if (length(A) == 1L) {
    A <- rep.int(as.integer(A), N)
  } else {
    A <- as.integer(A)
    if (length(A) != N) {
      stop("A must be NULL, length 1, or length N.", call. = FALSE)
    }
  }

  X1 <- stats::rnorm(N)
  X2 <- stats::rnorm(N)
  X3 <- stats::rbinom(N, size = 1L, prob = 0.5)
  X <- cbind(X1 = X1, X2 = X2, X3 = X3)

  scale <- get_prognosis_scale(prognosis, params)
  alpha <- get_alpha_effect(effect, params)
  heterogeneity_scale <- get_heterogeneity_scale(effect, effect_heterogeneity, params)
  linear_betas <- get_linear_endpoint_betas(endpoint_prognosis_pattern, params)

  zeta_D1 <- rep(params$theta_D10, N) +
    as.vector(X %*% (scale * linear_betas$D1))
  zeta_D2 <- rep(params$theta_D20, N) +
    as.vector(X %*% (scale * linear_betas$D2))
  zeta_D3 <- rep(params$theta_D30, N) +
    as.vector(X %*% (scale * linear_betas$D3))

  if (identical(outcome_shape, "nonlinear")) {
    nonlinear <- make_nonlinear_terms(X1, X3)
    zeta_D1 <- zeta_D1 + as.vector(nonlinear %*% (scale * params$beta_D1_nonlinear_base))
    zeta_D2 <- zeta_D2 + as.vector(nonlinear %*% (scale * params$beta_D2_nonlinear_base))
    zeta_D3 <- zeta_D3 + as.vector(nonlinear %*% (scale * params$beta_D3_nonlinear_base))
  }

  treatment_interaction_D1 <- heterogeneity_scale * as.vector((A * X) %*% params$gamma_D1_base)
  treatment_interaction_D2 <- heterogeneity_scale * as.vector((A * X) %*% params$gamma_D2_base)
  treatment_interaction_D3 <- heterogeneity_scale * as.vector((A * X) %*% params$gamma_D3_base)

  zeta_D1 <- zeta_D1 + alpha[["D1"]] * A + treatment_interaction_D1
  zeta_D2 <- zeta_D2 + alpha[["D2"]] * A + treatment_interaction_D2
  zeta_D3 <- zeta_D3 + alpha[["D3"]] * A + treatment_interaction_D3

  weibull_k <- params$weibull_k
  event_time <- (
    -weibull_k * log(stats::runif(N)) / (params$lambda0 * exp(zeta_D1))
  )^(1 / weibull_k)
  death_event <- as.integer(event_time <= params$S)
  death_time <- pmin(event_time, params$S)

  alsfrs_latent <- zeta_D2 + stats::rnorm(N, mean = 0, sd = params$sigma_D2)
  alsfrs <- pmin(48L, pmax(0L, as.integer(round(alsfrs_latent))))

  log_nfl <- zeta_D3 + stats::rnorm(N, mean = 0, sd = params$sigma_D3)
  nfl_score <- -log_nfl

  dat <- data.table::data.table(
    id = seq_len(N),
    A = A,
    X1 = X1,
    X2 = X2,
    X3 = X3,
    death_time = death_time,
    death_event = death_event,
    alsfrs = alsfrs,
    log_nfl = log_nfl,
    nfl_score = nfl_score,
    zeta_D1 = zeta_D1,
    zeta_D2 = zeta_D2,
    zeta_D3 = zeta_D3
  )

  if (n_noise_covariates >= 1L) {
    dat[, X4 := stats::rnorm(.N)]
  }
  if (n_noise_covariates >= 2L) {
    dat[, X5 := stats::rnorm(.N)]
  }
  adjustment_names <- if (identical(adjustment_covariate_pattern, "omit_X3_add_X4")) {
    c("X1", "X2", "X4")
  } else if (identical(adjustment_covariate_pattern, "omit_X3")) {
    c("X1", "X2")
  } else {
    intersect(paste0("X", seq_len(5L)), names(dat))
  }
  attr(dat, "adjustment_covariates") <- adjustment_names
  dat[]
}

generate_trial_data <- function(N,
                                e,
                                effect,
                                prognosis,
                                outcome_shape = "linear",
                                effect_heterogeneity = "homogeneous",
                                endpoint_prognosis_pattern = "shared",
                                n_noise_covariates = 0L,
                                adjustment_covariate_pattern = "all_available",
                                specification = NULL,
                                params = default_dgm_params()) {
  generate_subject_data(
    N = N,
    A = NULL,
    e = e,
    effect = effect,
    prognosis = prognosis,
    outcome_shape = outcome_shape,
    effect_heterogeneity = effect_heterogeneity,
    endpoint_prognosis_pattern = endpoint_prognosis_pattern,
    n_noise_covariates = n_noise_covariates,
    adjustment_covariate_pattern = adjustment_covariate_pattern,
    specification = specification,
    params = params
  )
}
