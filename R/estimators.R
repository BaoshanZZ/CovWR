clip_tau <- function(x, params = default_dgm_params()) {
  pmin(1 - params$min_tau, pmax(params$min_tau, x))
}

calc_measures <- function(tau_w, tau_l, params = default_dgm_params()) {
  tau_w <- clip_tau(tau_w, params)
  tau_l <- clip_tau(tau_l, params)
  delta <- tau_w - tau_l
  delta <- pmin(1 - params$min_tau, pmax(-1 + params$min_tau, delta))

  c(
    NB = tau_w - tau_l,
    LWR = log(tau_w) - log(tau_l),
    LWO = log1p(delta) - log1p(-delta)
  )
}

measure_gradient <- function(measure, tau_w, tau_l, params = default_dgm_params()) {
  measure <- match.arg(measure, c("NB", "LWR", "LWO"))
  tau_w <- clip_tau(tau_w, params)
  tau_l <- clip_tau(tau_l, params)
  delta <- tau_w - tau_l
  delta <- pmin(1 - params$min_tau, pmax(-1 + params$min_tau, delta))

  switch(
    measure,
    NB = c(w = 1, l = -1),
    LWR = c(w = 1 / tau_w, l = -1 / tau_l),
    LWO = c(w = 2 / (1 - delta^2), l = -2 / (1 - delta^2))
  )
}

wald_row <- function(method,
                     estimand,
                     estimate,
                     se,
                     truth,
                     null_value = 0,
                     alpha = 0.05) {
  zcrit <- stats::qnorm(1 - alpha / 2)
  ci_low <- estimate - zcrit * se
  ci_high <- estimate + zcrit * se
  z <- (estimate - null_value) / se
  p_value <- 2 * stats::pnorm(-abs(z))

  data.table::data.table(
    method = method,
    estimand = estimand,
    estimate = estimate,
    se = se,
    ci_low = ci_low,
    ci_high = ci_high,
    p_value = p_value,
    reject = p_value < alpha,
    truth = truth,
    covered = ci_low <= truth & truth <= ci_high
  )
}

if_se_for_measure <- function(psi_components, tau_w, tau_l, measure, N, params = default_dgm_params()) {
  grad <- measure_gradient(measure, tau_w, tau_l, params)
  psi_m <- as.vector(psi_components[, "w"] * grad[["w"]] + psi_components[, "l"] * grad[["l"]])
  sqrt(mean((psi_m - mean(psi_m))^2) / N)
}

fit_if_projection_components <- function(dat, pairwise, feature_fun = adjustment_features) {
  Fmat <- feature_fun(dat)
  Z <- dat$A - pairwise$ehat
  design <- Fmat * as.numeric(Z)

  fit_one <- function(y) {
    fit <- stats::lm.fit(x = design, y = y)
    beta <- fit$coefficients
    beta[is.na(beta)] <- 0
    beta
  }

  beta_w <- fit_one(pairwise$psi[, "w"])
  beta_l <- fit_one(pairwise$psi[, "l"])
  b_hat <- cbind(
    w = as.vector(Fmat %*% beta_w),
    l = as.vector(Fmat %*% beta_l)
  )
  psi_adj <- pairwise$psi - b_hat * as.numeric(Z)

  tau_adj <- c(
    w = pairwise$tau_w - mean(Z * b_hat[, "w"]),
    l = pairwise$tau_l - mean(Z * b_hat[, "l"])
  )

  list(
    beta = cbind(w = beta_w, l = beta_l),
    b_hat = b_hat,
    psi_adj = psi_adj,
    tau_adj = tau_adj
  )
}

fadj_feature_matrix <- function(dat,
                                learner = c("zero", "linear", "interactions", "splines")) {
  learner <- match.arg(learner)
  n <- nrow(dat)
  if (identical(learner, "zero")) {
    return(matrix(0, nrow = n, ncol = 0L))
  }

  X1 <- dat$X1
  X2 <- dat$X2
  X3 <- dat$X3
  if (identical(learner, "linear")) {
    return(cbind(intercept = 1, X1 = X1, X2 = X2, X3 = X3))
  }

  if (identical(learner, "interactions")) {
    return(cbind(
      intercept = 1,
      X1 = X1,
      X2 = X2,
      X3 = X3,
      X1_X2 = X1 * X2,
      X1_X3 = X1 * X3,
      X2_X3 = X2 * X3
    ))
  }

  sx1 <- splines::ns(X1, df = 3)
  sx2 <- splines::ns(X2, df = 3)
  colnames(sx1) <- paste0("ns_X1_", seq_len(ncol(sx1)))
  colnames(sx2) <- paste0("ns_X2_", seq_len(ncol(sx2)))
  cbind(intercept = 1, sx1, sx2, X3 = X3)
}

fit_direct_risk_linear <- function(Fmat, Z, y, train_idx = seq_along(y)) {
  if (ncol(Fmat) == 0L) {
    return(rep(0, ncol(Fmat)))
  }
  design <- Fmat[train_idx, , drop = FALSE] * as.numeric(Z[train_idx])
  fit <- stats::lm.fit(x = design, y = y[train_idx])
  beta <- fit$coefficients
  beta[is.na(beta)] <- 0
  beta
}

predict_direct_risk_linear <- function(Fmat, beta) {
  if (ncol(Fmat) == 0L) {
    return(rep(0, nrow(Fmat)))
  }
  as.vector(Fmat %*% beta)
}

make_cv_folds <- function(n, n_folds = 5L) {
  n_folds <- max(2L, min(as.integer(n_folds), n))
  ((seq_len(n) - 1L) %% n_folds) + 1L
}

select_fadj_learner <- function(dat,
                                y,
                                Z,
                                learners = c("zero", "linear", "interactions", "splines"),
                                n_folds = 5L) {
  folds <- make_cv_folds(nrow(dat), n_folds)
  risks <- vapply(learners, function(learner) {
    Fmat <- fadj_feature_matrix(dat, learner)
    fold_risk <- vapply(sort(unique(folds)), function(fold) {
      valid_idx <- which(folds == fold)
      train_idx <- which(folds != fold)
      beta <- fit_direct_risk_linear(Fmat, Z, y, train_idx)
      b_valid <- predict_direct_risk_linear(Fmat[valid_idx, , drop = FALSE], beta)
      mean((y[valid_idx] - Z[valid_idx] * b_valid)^2)
    }, numeric(1))
    mean(fold_risk)
  }, numeric(1))

  # Ties are resolved in learner order, so zero is selected unless another learner
  # strictly improves the direct-risk estimate.
  learner <- learners[which.min(risks)]
  list(learner = learner, risks = risks)
}

fit_fadj_components <- function(dat,
                                pairwise,
                                learners = c("zero", "linear", "interactions", "splines"),
                                n_folds = 5L) {
  Z <- dat$A - pairwise$ehat

  fit_one <- function(y) {
    selected <- select_fadj_learner(dat, y, Z, learners, n_folds)
    Fmat <- fadj_feature_matrix(dat, selected$learner)
    beta <- fit_direct_risk_linear(Fmat, Z, y)
    b_hat <- predict_direct_risk_linear(Fmat, beta)
    list(
      learner = selected$learner,
      cv_risks = selected$risks,
      beta = beta,
      b_hat = b_hat
    )
  }

  fit_w <- fit_one(pairwise$psi[, "w"])
  fit_l <- fit_one(pairwise$psi[, "l"])
  b_hat <- cbind(w = fit_w$b_hat, l = fit_l$b_hat)
  psi_adj <- pairwise$psi - b_hat * as.numeric(Z)

  tau_adj <- c(
    w = pairwise$tau_w - mean(Z * b_hat[, "w"]),
    l = pairwise$tau_l - mean(Z * b_hat[, "l"])
  )

  list(
    learners = c(w = fit_w$learner, l = fit_l$learner),
    cv_risks = list(w = fit_w$cv_risks, l = fit_l$cv_risks),
    b_hat = b_hat,
    psi_adj = psi_adj,
    tau_adj = tau_adj
  )
}

compute_sbw_weights <- function(X_arm, target_means, nonnegative = FALSE) {
  X_arm <- as.matrix(X_arm)
  n <- nrow(X_arm)
  if (n <= ncol(X_arm) + 1L) {
    stop("Not enough subjects in an arm to solve SBW constraints.", call. = FALSE)
  }

  eq_constraints <- cbind(intercept = rep(1, n), X_arm)
  eq_targets <- c(n, n * target_means)

  if (!nonnegative) {
    # Closed-form solution for exact equality SBW without non-negativity constraints:
    # minimize ||w - 1||^2 subject to C'w = target.
    lhs <- crossprod(eq_constraints)
    rhs <- eq_targets - colSums(eq_constraints)
    coef <- tryCatch(
      qr.solve(lhs, rhs),
      error = function(e) {
        ridge <- 1e-8 * diag(ncol(lhs))
        solve(lhs + ridge, rhs)
      }
    )
    return(as.numeric(1 + eq_constraints %*% coef))
  }

  if (nonnegative) {
    Amat <- cbind(eq_constraints, diag(n))
    bvec <- c(eq_targets, rep(0, n))
    meq <- ncol(eq_constraints)
  } else {
    Amat <- eq_constraints
    bvec <- eq_targets
    meq <- ncol(eq_constraints)
  }

  qp <- quadprog::solve.QP(
    Dmat = diag(n),
    dvec = rep(1, n),
    Amat = Amat,
    bvec = bvec,
    meq = meq
  )
  as.numeric(qp$solution)
}

estimate_sbw_tau <- function(dat, pairwise, nonnegative = FALSE) {
  X <- as.matrix(dat[, .(X1, X2, X3)])
  target <- colMeans(X)

  idx_t <- which(dat$A == 1L)
  idx_c <- which(dat$A == 0L)
  w_t <- compute_sbw_weights(X[idx_t, , drop = FALSE], target, nonnegative = nonnegative)
  w_c <- compute_sbw_weights(X[idx_c, , drop = FALSE], target, nonnegative = nonnegative)

  pl <- pairwise$placement
  delta_w <- weighted.mean(pl$g_w_T[idx_t], w_t) + weighted.mean(pl$g_w_C[idx_c], w_c)
  delta_l <- weighted.mean(pl$g_l_T[idx_t], w_t) + weighted.mean(pl$g_l_C[idx_c], w_c)

  list(
    tau = c(w = pairwise$tau_w + delta_w, l = pairwise$tau_l + delta_l),
    weights = list(treated = w_t, control = w_c)
  )
}

estimate_all_methods <- function(dat,
                                 truth_values = c(NB = NA_real_, LWR = NA_real_, LWO = NA_real_),
                                 params = default_dgm_params(),
                                 sbw_nonnegative = FALSE,
                                 include_sbw = FALSE,
                                 alpha = 0.05) {
  pairwise <- compute_pairwise_quantities(dat, params)
  N <- nrow(dat)
  rows <- list()

  add_method_rows <- function(method, tau, psi_components) {
    measures <- calc_measures(tau[["w"]], tau[["l"]], params)
    data.table::rbindlist(lapply(names(measures), function(m) {
      se <- if_se_for_measure(psi_components, tau[["w"]], tau[["l"]], m, N, params)
      wald_row(
        method = method,
        estimand = m,
        estimate = unname(measures[[m]]),
        se = se,
        truth = unname(truth_values[[m]]),
        null_value = 0,
        alpha = alpha
      )
    }))
  }

  tau_unadj <- c(w = pairwise$tau_w, l = pairwise$tau_l)
  rows[["UNADJ"]] <- add_method_rows("UNADJ", tau_unadj, pairwise$psi)

  ladj <- fit_if_projection_components(dat, pairwise)
  rows[["LADJ"]] <- add_method_rows("LADJ", ladj$tau_adj, ladj$psi_adj)

  fadj <- fit_fadj_components(dat, pairwise)
  rows[["FADJ"]] <- add_method_rows("FADJ", fadj$tau_adj, fadj$psi_adj)

  if (isTRUE(include_sbw)) {
    sbw <- estimate_sbw_tau(dat, pairwise, nonnegative = sbw_nonnegative)
    # The SBW standard error uses the first-order equivalent IF projection variance.
    rows[["SBW"]] <- add_method_rows("SBW", sbw$tau, ladj$psi_adj)
  }

  out <- data.table::rbindlist(rows, use.names = TRUE)
  tau_lookup <- list(
    UNADJ = tau_unadj,
    LADJ = ladj$tau_adj,
    FADJ = fadj$tau_adj
  )
  if (isTRUE(include_sbw)) {
    tau_lookup[["SBW"]] <- sbw$tau
  }
  out[, `:=`(
    tau_w = vapply(as.character(method), function(x) tau_lookup[[x]][["w"]], numeric(1)),
    tau_l = vapply(as.character(method), function(x) tau_lookup[[x]][["l"]], numeric(1)),
    learner_w = ifelse(method == "FADJ", unname(fadj$learners[["w"]]), NA_character_),
    learner_l = ifelse(method == "FADJ", unname(fadj$learners[["l"]]), NA_character_)
  )]
  out[]
}
