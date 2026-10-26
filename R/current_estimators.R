# Mechanically extracted from the validated retained implementation; estimator logic is unchanged.
fit_projection_qr <- function (X, A, psi, e, center = TRUE)
{
    X <- as.matrix(X)
    if (center)
        X <- sweep(X, 2, colMeans(X), "-")
    if (!ncol(X))
        return(list(beta = matrix(numeric(), 0, 2), adjustment = psi *
            0, residual = psi, rank = 0L, condition = 1, normal_error = 0))
    design <- X * (A - e)
    gram <- crossprod(design)
    qr_fit <- qr(design, tol = 1e-10)
    if (qr_fit$rank != ncol(X))
        stop("Rank-deficient LADJ design; no regularization applied.")
    rc <- rcond(gram)
    if (!is.finite(rc) || rc <= .Machine$double.eps)
        stop("Numerically singular LADJ Gram matrix.")
    beta <- qr.coef(qr_fit, psi)
    adjustment <- design %*% beta
    residual <- psi - adjustment
    list(beta = beta, adjustment = adjustment, residual = residual,
        rank = qr_fit$rank, condition = kappa(gram, exact = TRUE),
        normal_error = max(abs(crossprod(design, residual)))/nrow(X))
}

failed_estimator_rows <- function (W, reason)
{
    data.table::data.table(method = if (W == "W0")
        "UNADJ"
    else "LADJ", estimand = c("NB", "LWR", "LWO"), W = W, p = length(working_sets[[W]]),
        success = FALSE, failure_reason = reason, estimate = NA_real_,
        se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        p_value = NA_real_, reject = NA, truth = NA_real_, covered = NA,
        tau_w = NA_real_, tau_l = NA_real_, gram_rank = NA_integer_,
        gram_condition = NA_real_, normal_error = NA_real_)
}

run_win_estimators_for_all_W <- function (dat, params, truth, sets = working_sets)
{
    pair <- legacy_core$compute_pairwise_quantities(dat, params$e,
        params)
    tau <- c(w = pair$tau_w, l = pair$tau_l)
    rows <- lapply(names(sets), function(W) {
        tryCatch({
            cols <- names(population_means)[sets[[W]]]
            X <- if (length(cols))
                as.matrix(dat[, ..cols])
            else matrix(numeric(), nrow(dat), 0)
            fit <- fit_projection_qr(X, dat$A, pair$psi, params$e)
            tau_adj <- tau - colMeans(fit$adjustment)
            estimates <- legacy_core$calc_measures(tau_adj[[1]],
                tau_adj[[2]], params)
            out <- data.table::rbindlist(lapply(names(estimates),
                function(m) {
                  se <- legacy_core$if_se_for_measure(fit$residual,
                    tau_adj[[1]], tau_adj[[2]], m, nrow(dat),
                    params)
                  if (!is.finite(se) || se <= 0)
                    stop("Nonpositive or nonfinite standard error.")
                  legacy_core$wald_row(if (W == "W0")
                    "UNADJ"
                  else "LADJ", m, estimates[[m]], se, truth[[m]])
                }))
            out[, `:=`(W = W, p = length(cols), success = TRUE,
                failure_reason = "", tau_w = tau_adj[1], tau_l = tau_adj[2],
                gram_rank = fit$rank, gram_condition = fit$condition,
                normal_error = fit$normal_error, unadjusted_tau_w = tau[1],
                unadjusted_tau_l = tau[2])]
            out
        }, error = function(err) {
            failed <- failed_estimator_rows(W, conditionMessage(err))
            failed[, `:=`(unadjusted_tau_w = tau[1], unadjusted_tau_l = tau[2])]
            failed
        })
    })
    data.table::rbindlist(rows, fill = TRUE)
}
