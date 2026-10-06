validate_measure_domain <- function(tau_w, tau_l) {
  tau <- c(w = tau_w, l = tau_l)
  if (length(tau_w) != 1L || length(tau_l) != 1L || any(!is.finite(tau))) {
    stop("tau_w and tau_l must be finite scalar values.", call. = FALSE)
  }
  if (tau_w <= 0 || tau_l <= 0) {
    stop("LWR is undefined unless tau_w and tau_l are both positive.", call. = FALSE)
  }
  delta <- tau_w - tau_l
  if (abs(delta) >= 1) {
    stop("LWO is undefined unless abs(tau_w - tau_l) is less than one.", call. = FALSE)
  }
  invisible(tau)
}

calc_measures <- function(tau_w, tau_l, params = default_dgm_params()) {
  validate_measure_domain(tau_w, tau_l)
  delta <- tau_w - tau_l

  c(
    NB = tau_w - tau_l,
    LWR = log(tau_w) - log(tau_l),
    LWO = log1p(delta) - log1p(-delta)
  )
}

measure_gradient <- function(measure, tau_w, tau_l, params = default_dgm_params()) {
  measure <- match.arg(measure, c("NB", "LWR", "LWO"))
  validate_measure_domain(tau_w, tau_l)
  delta <- tau_w - tau_l

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

add_bootstrap_wald_columns <- function(out,
                                       alpha = 0.05,
                                       null_value = 0) {
  zcrit <- stats::qnorm(1 - alpha / 2)
  out[, `:=`(
    ci_low_boot = estimate - zcrit * se_boot,
    ci_high_boot = estimate + zcrit * se_boot,
    p_value_boot = 2 * stats::pnorm(-abs((estimate - null_value) / se_boot))
  )]
  out[, `:=`(
    reject_boot = p_value_boot < alpha,
    covered_boot = ci_low_boot <= truth & truth <= ci_high_boot
  )]
  out[]
}

if_se_for_measure <- function(psi_components, tau_w, tau_l, measure, N, params = default_dgm_params()) {
  grad <- measure_gradient(measure, tau_w, tau_l, params)
  psi_mat <- as.matrix(psi_components[, c("w", "l"), drop = FALSE])
  centered <- sweep(psi_mat, 2L, colMeans(psi_mat), FUN = "-")
  sigma_hat <- crossprod(centered) / N
  grad_vec <- matrix(c(grad[["w"]], grad[["l"]]), ncol = 1L)
  var_hat <- as.numeric(crossprod(grad_vec, sigma_hat %*% grad_vec) / N)
  sqrt(max(var_hat, 0))
}

fit_if_projection_components <- function(dat, pairwise) {
  Fmat <- adjustment_features(dat)
  Z <- dat$A - pairwise$allocation_prob
  design <- Fmat * as.numeric(Z)
  gram <- crossprod(design)
  gram_rank <- qr(gram)$rank
  if (gram_rank < ncol(gram)) {
    stop(
      "The LADJ Gram matrix is singular (rank ", gram_rank,
      " of ", ncol(gram), ").",
      call. = FALSE
    )
  }
  reciprocal_condition <- rcond(gram)
  if (!is.finite(reciprocal_condition) || reciprocal_condition <= .Machine$double.eps) {
    stop("The LADJ Gram matrix is numerically singular.", call. = FALSE)
  }

  fit_one <- function(y) {
    score <- crossprod(design, y)
    beta <- solve(gram, score)
    stats::setNames(as.numeric(beta), colnames(Fmat))
  }

  beta_w <- fit_one(pairwise$psi[, "w"])
  beta_l <- fit_one(pairwise$psi[, "l"])
  b_hat <- cbind(
    w = as.vector(Fmat %*% beta_w),
    l = as.vector(Fmat %*% beta_l)
  )
  psi_adj <- pairwise$psi - b_hat * as.numeric(Z)

  idx_t <- dat$A == 1L
  idx_c <- !idx_t
  imbalance <- colMeans(Fmat[idx_t, , drop = FALSE]) -
    colMeans(Fmat[idx_c, , drop = FALSE])
  imbalance_factor <- sum(idx_t) * sum(idx_c) / nrow(dat)^2
  tau_adj <- c(
    w = pairwise$tau_w - imbalance_factor * sum(beta_w * imbalance),
    l = pairwise$tau_l - imbalance_factor * sum(beta_l * imbalance)
  )

  list(
    beta = cbind(w = beta_w, l = beta_l),
    b_hat = b_hat,
    psi_adj = psi_adj,
    tau_adj = tau_adj,
    imbalance = imbalance,
    imbalance_factor = imbalance_factor,
    gram_rank = gram_rank,
    reciprocal_condition = reciprocal_condition
  )
}

fadj_feature_matrix <- function(dat,
                                learner = c("zero", "linear", "interactions", "splines")) {
  learner <- match.arg(learner)
  n <- nrow(dat)
  if (identical(learner, "zero")) {
    return(matrix(0, nrow = n, ncol = 0L))
  }

  covariate_names <- adjustment_covariate_names(dat)
  X <- as.matrix(dat[, ..covariate_names])
  if (identical(learner, "linear")) {
    return(cbind(intercept = 1, X))
  }

  if (identical(learner, "interactions")) {
    pairs <- utils::combn(covariate_names, 2L)
    interaction_terms <- vapply(seq_len(ncol(pairs)), function(j) {
      X[, pairs[1L, j]] * X[, pairs[2L, j]]
    }, numeric(n))
    colnames(interaction_terms) <- apply(pairs, 2L, paste, collapse = "_")
    return(cbind(intercept = 1, X, interaction_terms))
  }

  continuous_names <- intersect(c("X1", "X2", "X4", "X5"), covariate_names)
  spline_terms <- lapply(continuous_names, function(name) {
    out <- splines::ns(dat[[name]], df = 3)
    colnames(out) <- paste0("ns_", name, "_", seq_len(ncol(out)))
    out
  })
  binary_names <- intersect("X3", covariate_names)
  do.call(cbind, c(list(intercept = rep(1, n)), spline_terms, list(X[, binary_names, drop = FALSE])))
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
  Z <- dat$A - pairwise$allocation_prob

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

stratified_bootstrap_sample <- function(dat) {
  idx_t <- which(dat$A == 1L)
  idx_c <- which(dat$A == 0L)
  boot_idx <- c(
    sample(idx_t, length(idx_t), replace = TRUE),
    sample(idx_c, length(idx_c), replace = TRUE)
  )
  boot <- data.table::copy(dat[boot_idx])
  boot[, id := seq_len(.N)]
  boot[]
}

bootstrap_method_estimates <- function(dat,
                                       B,
                                       allocation_prob,
                                       params = default_dgm_params(),
                                       include_sbw = FALSE,
                                       sbw_nonnegative = FALSE,
                                       method_names = c("UNADJ", "LADJ"),
                                       alpha = 0.05) {
  if (B <= 0L) {
    return(NULL)
  }

  truth_na <- c(NB = NA_real_, LWR = NA_real_, LWO = NA_real_)
  rows <- vector("list", B)
  for (b in seq_len(B)) {
    boot_dat <- stratified_bootstrap_sample(dat)
    boot_est <- estimate_all_methods(
      dat = boot_dat,
      allocation_prob = allocation_prob,
      truth_values = truth_na,
      params = params,
      sbw_nonnegative = sbw_nonnegative,
      include_sbw = include_sbw,
      bootstrap_B = 0L,
      compute_inference = FALSE,
      method_names = method_names,
      alpha = alpha
    )
    rows[[b]] <- boot_est[, .(
      bootstrap_id = b,
      method,
      estimand,
      estimate
    )]
  }

  boot <- data.table::rbindlist(rows, use.names = TRUE)
  boot[, .(
    bootstrap_B = .N,
    boot_var = stats::var(estimate),
    se_boot = stats::sd(estimate)
  ), by = .(method, estimand)]
}

estimate_all_methods <- function(dat,
                                 allocation_prob,
                                 truth_values = c(NB = NA_real_, LWR = NA_real_, LWO = NA_real_),
                                 params = default_dgm_params(),
                                 sbw_nonnegative = FALSE,
                                 include_sbw = FALSE,
                                 bootstrap_B = 0L,
                                 compute_inference = TRUE,
                                 method_names = c("UNADJ", "LADJ"),
                                 alpha = 0.05) {
  method_names <- unique(toupper(method_names))
  allowed_methods <- c("UNADJ", "LADJ", "FADJ", "SBW")
  bad_methods <- setdiff(method_names, allowed_methods)
  if (length(bad_methods) > 0L) {
    stop("Unknown method_names value(s): ", paste(bad_methods, collapse = ", "), call. = FALSE)
  }
  if ("SBW" %in% method_names) {
    include_sbw <- TRUE
  }

  pairwise <- compute_pairwise_quantities(
    dat = dat,
    allocation_prob = allocation_prob,
    params = params
  )
  N <- nrow(dat)
  rows <- list()
  tau_lookup <- list()
  fadj <- NULL

  add_method_rows <- function(method, tau, psi_components) {
    measures <- calc_measures(tau[["w"]], tau[["l"]], params)
    data.table::rbindlist(lapply(names(measures), function(m) {
      if (!isTRUE(compute_inference)) {
        return(data.table::data.table(
          method = method,
          estimand = m,
          estimate = unname(measures[[m]])
        ))
      }
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

  if ("UNADJ" %in% method_names) {
    tau_unadj <- c(w = pairwise$tau_w, l = pairwise$tau_l)
    rows[["UNADJ"]] <- add_method_rows("UNADJ", tau_unadj, pairwise$psi)
    tau_lookup[["UNADJ"]] <- tau_unadj
  }

  needs_ladj <- any(c("LADJ", "SBW") %in% method_names)
  if (needs_ladj) {
    ladj <- fit_if_projection_components(dat, pairwise)
    if ("LADJ" %in% method_names) {
      rows[["LADJ"]] <- add_method_rows("LADJ", ladj$tau_adj, ladj$psi_adj)
      tau_lookup[["LADJ"]] <- ladj$tau_adj
    }
  }

  if ("FADJ" %in% method_names) {
    fadj <- fit_fadj_components(dat, pairwise)
    rows[["FADJ"]] <- add_method_rows("FADJ", fadj$tau_adj, fadj$psi_adj)
    tau_lookup[["FADJ"]] <- fadj$tau_adj
  }

  if (isTRUE(include_sbw) && "SBW" %in% method_names) {
    sbw <- estimate_sbw_tau(dat, pairwise, nonnegative = sbw_nonnegative)
    # The SBW standard error uses the first-order equivalent IF projection variance.
    rows[["SBW"]] <- add_method_rows("SBW", sbw$tau, ladj$psi_adj)
    tau_lookup[["SBW"]] <- sbw$tau
  }

  out <- data.table::rbindlist(rows, use.names = TRUE)
  out[, `:=`(
    tau_w = vapply(as.character(method), function(x) tau_lookup[[x]][["w"]], numeric(1)),
    tau_l = vapply(as.character(method), function(x) tau_lookup[[x]][["l"]], numeric(1)),
    learner_w = if (!is.null(fadj)) ifelse(method == "FADJ", unname(fadj$learners[["w"]]), NA_character_) else NA_character_,
    learner_l = if (!is.null(fadj)) ifelse(method == "FADJ", unname(fadj$learners[["l"]]), NA_character_) else NA_character_
  )]

  if (bootstrap_B > 0L) {
    boot_se <- bootstrap_method_estimates(
      dat = dat,
      B = as.integer(bootstrap_B),
      allocation_prob = allocation_prob,
      params = params,
      include_sbw = include_sbw,
      sbw_nonnegative = sbw_nonnegative,
      method_names = method_names,
      alpha = alpha
    )
    out <- merge(out, boot_se, by = c("method", "estimand"), all.x = TRUE, sort = FALSE)
    out <- add_bootstrap_wald_columns(out, alpha = alpha)
  } else {
    out[, `:=`(
      bootstrap_B = as.integer(bootstrap_B),
      boot_var = NA_real_,
      se_boot = NA_real_,
      ci_low_boot = NA_real_,
      ci_high_boot = NA_real_,
      p_value_boot = NA_real_,
      reject_boot = NA,
      covered_boot = NA
    )]
  }

  out[]
}
