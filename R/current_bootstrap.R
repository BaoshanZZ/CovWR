# Mechanically extracted from the validated retained implementation; estimator logic is unchanged.
bootstrap_context <- function (dat, p)
{
    pair <- legacy_core$compute_pairwise_quantities(dat, p$e,
        p)
    list(dat = dat, X = as.matrix(dat[, names(population_means),
        with = FALSE]), A = dat$A, N = nrow(dat), ti = which(dat$A ==
        1), ci = which(dat$A == 0), win = pair$matrices$win,
        loss = pair$matrices$loss, p = p)
}

bootstrap_placements <- function (ctx, counts)
{
    if (is.null(dim(counts)))
        counts <- matrix(counts, ncol = 1)
    ct <- counts[ctx$ti, , drop = FALSE]
    cc <- counts[ctx$ci, , drop = FALSE]
    n1 <- colSums(ct)
    n0 <- colSums(cc)
    wt <- sweep(ctx$win %*% cc, 2, pmax(n0, 1), "/")
    lt <- sweep(ctx$loss %*% cc, 2, pmax(n0, 1), "/")
    wc <- sweep(t(ctx$win) %*% ct, 2, pmax(n1, 1), "/")
    lc <- sweep(t(ctx$loss) %*% ct, 2, pmax(n1, 1), "/")
    list(wt = wt, lt = lt, wc = wc, lc = lc, n1 = n1, n0 = n0,
        tau_w = colSums(ct * wt)/pmax(n1, 1), tau_l = colSums(ct *
            lt)/pmax(n1, 1))
}

bootstrap_refit <- function (ctx, counts, pl, b)
{
    sets <- working_sets
    tau <- c(pl$tau_w[b], pl$tau_l[b])
    estimates <- matrix(NA_real_, length(sets), 3, dimnames = list(names(sets),
        c("NB", "LWR", "LWO")))
    taus <- matrix(NA_real_, length(sets), 2)
    reason <- detail <- rep("", length(sets))
    ranks <- rep(NA_integer_, length(sets))
    conditions <- rep(NA_real_, length(sets))
    if (pl$n1[b] == 0 || pl$n0[b] == 0) {
        reason[] <- "missing_treatment_arm"
        return(list(estimates = estimates, tau = taus, reason = reason,
            detail = detail, rank = ranks, condition = conditions))
    }
    psi <- matrix(0, ctx$N, 2)
    psi[ctx$ti, ] <- cbind(pl$wt[, b] - tau[1], pl$lt[, b] -
        tau[2])/ctx$p$e
    psi[ctx$ci, ] <- cbind(pl$wc[, b] - tau[1], pl$lc[, b] -
        tau[2])/(1 - ctx$p$e)
    keep <- counts > 0
    Xc <- sweep(ctx$X, 2, as.vector(crossprod(counts, ctx$X))/ctx$N,
        "-")
    design <- Xc * (ctx$A - ctx$p$e)
    imbalance <- as.vector(crossprod(counts, design))/ctx$N
    root_counts <- sqrt(counts[keep])
    weighted_design <- design[keep, , drop = FALSE] * root_counts
    weighted_psi <- psi[keep, , drop = FALSE] * root_counts
    for (w in seq_along(sets)) {
        cols <- sets[[w]]
        ta <- tau
        if (length(cols)) {
            D <- weighted_design[, cols, drop = FALSE]
            qr_fit <- qr(D, tol = 1e-10)
            ranks[w] <- qr_fit$rank
            if (qr_fit$rank < length(cols)) {
                reason[w] <- "rank_deficient_projection"
                detail[w] <- paste("rank", qr_fit$rank, "p",
                  length(cols))
                next
            }
            gram <- crossprod(D)
            rc <- rcond(gram)
            if (!is.finite(rc) || rc <= .Machine$double.eps) {
                reason[w] <- "numerically_singular_gram"
                detail[w] <- paste("rcond", rc)
                next
            }
            conditions[w] <- kappa(gram, exact = TRUE)
            beta <- qr.coef(qr_fit, weighted_psi)
            ta <- tau - as.vector(crossprod(imbalance[cols],
                beta))
        }
        else {
            ranks[w] <- 0L
            conditions[w] <- 1
        }
        taus[w, ] <- ta
        value <- tryCatch(legacy_core$calc_measures(ta[1], ta[2],
            ctx$p), error = function(e) e)
        if (inherits(value, "error")) {
            reason[w] <- "invalid_measure_domain"
            detail[w] <- conditionMessage(value)
        }
        else estimates[w, ] <- value
    }
    list(estimates = estimates, tau = taus, reason = reason,
        detail = detail, rank = ranks, condition = conditions)
}
