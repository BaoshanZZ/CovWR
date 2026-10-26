# Mechanically extracted from the validated retained implementation; estimator logic is unchanged.
validate_measure_domain <- function (tau_w, tau_l)
{
    tau <- c(w = tau_w, l = tau_l)
    if (length(tau_w) != 1L || length(tau_l) != 1L || any(!is.finite(tau))) {
        stop("tau_w and tau_l must be finite scalar values.",
            call. = FALSE)
    }
    if (tau_w <= 0 || tau_l <= 0) {
        stop("LWR is undefined unless tau_w and tau_l are both positive.",
            call. = FALSE)
    }
    delta <- tau_w - tau_l
    if (abs(delta) >= 1) {
        stop("LWO is undefined unless abs(tau_w - tau_l) is less than one.",
            call. = FALSE)
    }
    invisible(tau)
}

calc_measures <- function (tau_w, tau_l, params = default_dgm_params())
{
    validate_measure_domain(tau_w, tau_l)
    delta <- tau_w - tau_l
    c(NB = tau_w - tau_l, LWR = log(tau_w) - log(tau_l), LWO = log1p(delta) -
        log1p(-delta))
}

measure_gradient <- function (measure, tau_w, tau_l, params = default_dgm_params())
{
    measure <- match.arg(measure, c("NB", "LWR", "LWO"))
    validate_measure_domain(tau_w, tau_l)
    delta <- tau_w - tau_l
    switch(measure, NB = c(w = 1, l = -1), LWR = c(w = 1/tau_w,
        l = -1/tau_l), LWO = c(w = 2/(1 - delta^2), l = -2/(1 -
        delta^2)))
}

wald_row <- function (method, estimand, estimate, se, truth, null_value = 0,
    alpha = 0.05)
{
    zcrit <- stats::qnorm(1 - alpha/2)
    ci_low <- estimate - zcrit * se
    ci_high <- estimate + zcrit * se
    z <- (estimate - null_value)/se
    p_value <- 2 * stats::pnorm(-abs(z))
    data.table::data.table(method = method, estimand = estimand,
        estimate = estimate, se = se, ci_low = ci_low, ci_high = ci_high,
        p_value = p_value, reject = p_value < alpha, truth = truth,
        covered = ci_low <= truth & truth <= ci_high)
}

if_se_for_measure <- function (psi_components, tau_w, tau_l, measure, N, params = default_dgm_params())
{
    grad <- measure_gradient(measure, tau_w, tau_l, params)
    psi_mat <- as.matrix(psi_components[, c("w", "l"), drop = FALSE])
    centered <- sweep(psi_mat, 2L, colMeans(psi_mat), FUN = "-")
    sigma_hat <- crossprod(centered)/N
    grad_vec <- matrix(c(grad[["w"]], grad[["l"]]), ncol = 1L)
    var_hat <- as.numeric(crossprod(grad_vec, sigma_hat %*% grad_vec)/N)
    sqrt(max(var_hat, 0))
}

compute_pairwise_quantities <- function (dat, allocation_prob, params = default_dgm_params())
{
    allocation_prob <- as.numeric(allocation_prob)
    if (length(allocation_prob) != 1L || is.na(allocation_prob) ||
        allocation_prob <= 0 || allocation_prob >= 1) {
        stop("allocation_prob must be a single known probability in (0, 1).",
            call. = FALSE)
    }
    treated <- dat[A == 1L]
    control <- dat[A == 0L]
    mats <- hierarchical_pair_matrices(treated, control, params)
    tau_w <- mean(mats$win)
    tau_l <- mean(mats$loss)
    g_w_T <- rowMeans(mats$win) - tau_w
    g_w_C <- colMeans(mats$win) - tau_w
    g_l_T <- rowMeans(mats$loss) - tau_l
    g_l_C <- colMeans(mats$loss) - tau_l
    placement <- data.table::data.table(id = c(treated$id, control$id),
        g_w_T = c(g_w_T, rep(NA_real_, nrow(control))), g_w_C = c(rep(NA_real_,
            nrow(treated)), g_w_C), g_l_T = c(g_l_T, rep(NA_real_,
            nrow(control))), g_l_C = c(rep(NA_real_, nrow(treated)),
            g_l_C))
    data.table::setkey(placement, id)
    placement <- placement[dat[, .(id)]]
    zero_na <- function(x) {
        x[is.na(x)] <- 0
        x
    }
    psi_w <- dat$A/allocation_prob * zero_na(placement$g_w_T) +
        (1 - dat$A)/(1 - allocation_prob) * zero_na(placement$g_w_C)
    psi_l <- dat$A/allocation_prob * zero_na(placement$g_l_T) +
        (1 - dat$A)/(1 - allocation_prob) * zero_na(placement$g_l_C)
    list(tau_w = tau_w, tau_l = tau_l, allocation_prob = allocation_prob,
        realized_allocation = mean(dat$A), placement = placement,
        psi = cbind(w = psi_w, l = psi_l), matrices = mats)
}

check_kernel_partition <- function (pair_obj)
{
    sums <- pair_obj$win + pair_obj$loss + pair_obj$tie
    isTRUE(all.equal(as.numeric(sums), rep(1, length(sums)),
        tolerance = 0))
}
