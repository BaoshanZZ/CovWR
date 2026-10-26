# Mechanically extracted from the validated retained implementation; estimator logic is unchanged.
full_population_draws <- function (n, A = 0L, e = 0.5)
{
    V <- sqrt(0.2) * rnorm(n) + sqrt(0.8) * matrix(rnorm(n *
        5L), n, 5L)
    X <- matrix(0, n, 10L, dimnames = list(NULL, names(population_means)[1:10]))
    X[, 1:5] <- V
    for (j in c(2L, 4L)) X[, j] <- as.numeric(V[, j] > qnorm(1 -
        population_means[j]))
    X[, c(6L, 8L, 10L)] <- matrix(rnorm(n * 3L), n, 3L)
    X[, 7L] <- rbinom(n, 1L, 0.5)
    X[, 9L] <- rbinom(n, 1L, 0.3)
    if (is.null(A))
        A <- rbinom(n, 1L, e)
    draws <- list(X = X, A = rep_len(A, n), U1 = runif(n), eps2 = rnorm(n),
        eps3 = rnorm(n))
    noise <- matrix(rnorm(n * 5L), n, 5L, dimnames = list(NULL,
        paste0("X", 11:15)))
    draws$X <- cbind(draws$X, noise)
    draws
}

staggered_predictors <- function (draws, p, effect = "null")
{
    beta <- cbind(p$pattern1, p$sigma2 * p$pattern2, p$sigma3 *
        p$pattern3)
    alpha <- if (effect == "null")
        rep(0, 3)
    else p$treatment_strength * p$alpha_base
    sweep(draws$X %*% beta, 2, c(p$theta1, p$theta2, p$theta3),
        "+") + outer(draws$A, alpha)
}

staggered_outcomes <- function (draws, p, effect = "null")
{
    zeta <- staggered_predictors(draws, p, effect)
    T1 <- (-p$shape * log(draws$U1)/(p$lambda0 * exp(zeta[, 1])))^(1/p$shape)
    event <- as.integer(T1 <= p$S)
    time <- pmin(T1, p$S)
    d2 <- pmin(p$S, pmax(0, round(zeta[, 2] + p$sigma2 * draws$eps2)))
    d2[event == 1] <- pmin(d2[event == 1], floor(time[event ==
        1]))
    log_d3 <- zeta[, 3] + p$sigma3 * draws$eps3
    log_d3[event == 1] <- NA_real_
    dat <- data.table::as.data.table(draws$X)
    dat[, `:=`(id = seq_len(.N), A = draws$A, death_time = time,
        death_event = event, icu_days = d2, log_biomarker = log_d3,
        biomarker = exp(log_d3))]
    dat
}
