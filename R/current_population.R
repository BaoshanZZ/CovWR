# Mechanically extracted from the validated retained implementation; estimator logic is unchanged.
binary_auc <- function (event, risk)
{
    n1 <- sum(event == 1)
    n0 <- sum(event == 0)
    if (!n1 || !n0)
        stop("AUC needs both event states.")
    (sum(rank(risk, ties.method = "average")[event == 1]) - n1 *
        (n1 + 1)/2)/(n1 * n0)
}

prediction_r2 <- function (y, pred)
1 - mean((y - pred)^2)/var(y)

survival_metrics <- function (dat, risk)
{
    obs <- data.frame(time = dat$death_time, event = dat$death_event,
        risk = risk)
    C <- survival::concordance(survival::Surv(time, event) ~
        risk, data = obs, reverse = TRUE, timewt = "n")$concordance
    c(Harrell_C = C, AUC_S = binary_auc(dat$death_event, risk))
}

expected_icu <- function (mu, sigma, S)
{
    ans <- numeric(length(mu))
    for (d in seq_len(S)) ans <- ans + pnorm((mu - (d - 0.5))/sigma)
    ans
}

calculate_population_win_truth <- function (params, settings, output_dir)
{
    streams <- make_streams(settings$seed_truth, 2)
    results <- lapply(seq_along(streams), function(k) {
        effect <- c("null", "alternative")[k]
        use_stream(streams[[k]])
        counts <- c(win = 0, loss = 0, tie = 0, D1 = 0, D2 = 0,
            D3 = 0)
        for (start in seq.int(1L, settings$truth_pairs, by = settings$truth_chunk)) {
            n <- min(settings$truth_chunk, settings$truth_pairs -
                start + 1L)
            t <- generate_hce_outcomes(make_population_draws(n,
                1L), params, effect)
            c <- generate_hce_outcomes(make_population_draws(n,
                0L), params, effect)
            pairs <- compare_hce_pair(t, c, params)
            counts[1:3] <- counts[1:3] + colSums(pairs[, .(win,
                loss, tie)])
            counts[4:6] <- counts[4:6] + as.numeric(table(factor(pairs$deciding_component,
                levels = c("D1", "D2", "D3"))))
        }
        prob <- counts/settings$truth_pairs
        tau <- unname(prob[1:2])
        raw_nb <- tau[1] - tau[2]
        if (effect == "null")
            tau[] <- mean(tau)
        measures <- legacy_core$calc_measures(tau[[1]], tau[[2]],
            params)
        sigma <- diag(unname(prob[1:2])) - tcrossprod(unname(prob[1:2]))
        ses <- vapply(c("NB", "LWR", "LWO"), function(m) {
            g <- legacy_core$measure_gradient(m, tau[[1]], tau[[2]],
                params)
            sqrt(as.numeric(crossprod(g, sigma %*% g))/settings$truth_pairs)
        }, numeric(1))
        data.table::data.table(effect = effect, n_pairs = settings$truth_pairs,
            tau_w = tau[1], tau_l = tau[2], tau_tie = prob[["tie"]],
            raw_tau_w = prob[["win"]], raw_tau_l = prob[["loss"]],
            raw_NB = raw_nb, raw_NB_mcse = ses[["NB"]], NB = measures[["NB"]],
            LWR = measures[["LWR"]], LWO = measures[["LWO"]],
            NB_truth_mcse = if (effect == "null")
                0
            else ses[["NB"]], LWR_truth_mcse = if (effect ==
                "null")
                0
            else ses[["LWR"]], LWO_truth_mcse = if (effect ==
                "null")
                0
            else ses[["LWO"]], resolved_D1 = prob[["D1"]], resolved_D2 = prob[["D2"]],
            resolved_D3 = prob[["D3"]], overall_tie = prob[["tie"]])
    })
    out <- data.table::rbindlist(results)
    save_csv(out, "hierarchy_calibration", output_dir)
    out
}

population_if_sample <- function (effect, params, n, reference_n, block)
{
    dat <- generate_rct_replication(n, effect, params)
    ref1 <- generate_hce_outcomes(make_population_draws(reference_n,
        1L), params, effect)
    ref0 <- generate_hce_outcomes(make_population_draws(reference_n,
        0L), params, effect)
    placements <- matrix(NA_real_, n, 2, dimnames = list(NULL,
        c("w", "l")))
    ti <- which(dat$A == 1)
    ci <- which(dat$A == 0)
    placements[ti, ] <- population_placements(dat[ti], ref0,
        params, block)
    placements[ci, ] <- population_placements(dat[ci], ref1,
        params, block)[, c("l", "w")]
    list(dat = dat, placements = placements)
}

calculate_oracle_if_pvr <- function (params, truth, settings, output_dir)
{
    streams <- make_streams(settings$seed_oracle, 4)
    results <- list()
    for (k in 1:2) {
        effect <- c("null", "alternative")[k]
        effect_value <- effect
        tau <- unlist(truth[effect == effect_value, .(tau_w,
            tau_l)], use.names = FALSE)
        use_stream(streams[[2 * k - 1L]])
        tr <- population_if_sample(effect, params, settings$oracle_train_n,
            settings$reference_n, settings$placement_block)
        use_stream(streams[[2 * k]])
        va <- population_if_sample(effect, params, settings$oracle_validation_n,
            settings$reference_n, settings$placement_block)
        score <- function(obj) {
            sweep(obj$placements, 2, tau, "-")/ifelse(obj$dat$A ==
                1, params$e, 1 - params$e)
        }
        psi_t <- score(tr)
        psi_v <- score(va)
        for (W in names(working_sets)) {
            cols <- names(population_means)[working_sets[[W]]]
            Xt <- if (length(cols))
                sweep(as.matrix(tr$dat[, ..cols]), 2, population_means[cols],
                  "-")
            else matrix(numeric(), nrow(tr$dat), 0)
            Xv <- if (length(cols))
                sweep(as.matrix(va$dat[, ..cols]), 2, population_means[cols],
                  "-")
            else matrix(numeric(), nrow(va$dat), 0)
            fit <- fit_projection_qr(Xt, tr$dat$A, psi_t, params$e,
                center = FALSE)
            residual <- psi_v - (Xv * (va$dat$A - params$e)) %*%
                fit$beta
            sigma0 <- cov(psi_v)
            sigma1 <- cov(residual)
            for (m in c("NB", "LWR", "LWO")) {
                g <- legacy_core$measure_gradient(m, tau[[1]],
                  tau[[2]], params)
                v0 <- as.numeric(crossprod(g, sigma0 %*% g))
                v1 <- as.numeric(crossprod(g, sigma1 %*% g))
                u <- as.vector(psi_v %*% g)
                r <- as.vector(residual %*% g)
                influence_ratio <- ((r - mean(r))^2 - (v1/v0) *
                  (u - mean(u))^2)/v0
                results[[length(results) + 1L]] <- data.table::data.table(effect = effect,
                  W = W, p = length(cols), estimand = m, V_unadjusted = v0,
                  V_adjusted = v1, oracle_PVR = 100 * (1 - v1/v0),
                  conditional_mcse_PVR = 100 * sd(influence_ratio)/sqrt(length(r)),
                  training_n = nrow(tr$dat), evaluation_n = nrow(va$dat),
                  reference_per_arm_per_split = settings$reference_n,
                  gram_rank = fit$rank, gram_condition = fit$condition)
            }
        }
        message("Population IF benchmark complete: ", effect)
    }
    out <- data.table::rbindlist(results)
    save_csv(out, "oracle_if_pvr_by_W", output_dir)
    out
}
