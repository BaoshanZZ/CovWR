# Mechanically extracted from the validated retained implementation; estimator logic is unchanged.
set_stream_seed <- function (seed)
{
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
}

make_streams <- function (seed, n)
{
    set_stream_seed(seed)
    out <- vector("list", n)
    s <- .Random.seed
    for (i in seq_len(n)) {
        out[[i]] <- s
        s <- parallel::nextRNGStream(s)
    }
    out
}

use_stream <- function (s)
assign(".Random.seed", s, envir = .GlobalEnv)

save_csv <- function (x, name, output_dir)
{
    data.table::fwrite(x, file.path(output_dir, paste0(name,
        ".csv")), na = "NA")
}

wilson_ci <- function (p, n)
{
    z <- qnorm(0.975)
    center <- (p + z^2/(2 * n))/(1 + z^2/n)
    half <- z * sqrt(p * (1 - p)/n + z^2/(4 * n^2))/(1 + z^2/n)
    cbind(low = center - half, high = center + half)
}

run_trial_tasks <- function (tasks, streams, params, truth, cores = 1L)
{
    run <- function(i) {
        task <- tasks[i]
        use_stream(streams[[i]])
        started <- proc.time()[[3]]
        tryCatch({
            tr <- truth[effect == task$effect]
            truths <- c(NB = tr$NB, LWR = tr$LWR, LWO = tr$LWO)
            dat <- generate_rct_replication(task$N, task$effect,
                params)
            raw <- run_win_estimators_for_all_W(dat, params,
                truths)
            raw[, `:=`(N = task$N, effect = task$effect, rep = task$rep,
                task_id = task$task_id, elapsed_seconds = proc.time()[[3]] -
                  started)]
            raw
        }, error = function(e) {
            out <- data.table::rbindlist(lapply(names(working_sets),
                failed_estimator_rows, reason = conditionMessage(e)))
            out[, `:=`(N = task$N, effect = task$effect, rep = task$rep,
                task_id = task$task_id, unadjusted_tau_w = NA_real_,
                unadjusted_tau_l = NA_real_, elapsed_seconds = proc.time()[[3]] -
                  started)]
            out
        })
    }
    ans <- if (cores > 1L)
        parallel::mclapply(seq_len(nrow(tasks)), run, mc.cores = cores,
            mc.set.seed = FALSE)
    else lapply(seq_len(nrow(tasks)), run)
    if (any(vapply(ans, inherits, logical(1), "try-error")))
        stop("Worker failed outside error capture.")
    data.table::rbindlist(ans, fill = TRUE)
}

summarize_operating_characteristics <- function (raw, oracle, output_dir)
{
    out <- raw[, {
        key <- .BY
        s <- .SD[success == TRUE]
        b <- raw[N == key$N & effect == key$effect & W == "W0" &
            estimand == key$estimand & success == TRUE]
        common <- merge(s[, .(rep, adj = estimate, se_adj = se)],
            b[, .(rep, unadj = estimate, se_unadj = se)], by = "rep")
        vu <- var(common$unadj)
        va <- var(common$adj)
        list(n_attempted = .N, n_success = nrow(s), n_failed = sum(!success),
            n_paired_with_unadjusted = nrow(common), estimate_mean = mean(s$estimate),
            truth = if (nrow(s)) s$truth[1] else NA_real_, bias = mean(s$estimate -
                s$truth), ESD = sd(s$estimate), ASE = mean(s$se),
            RMSSE = sqrt(mean(s$se^2)), ASE_ESD = mean(s$se)/sd(s$estimate),
            RMSSE_ESD = sqrt(mean(s$se^2)/var(s$estimate)), empirical_variance = var(s$estimate),
            mean_estimated_variance = mean(s$se^2), empirical_PVR = 100 *
                (1 - va/vu), estimated_PVR = 100 * (1 - mean(common$se_adj^2)/mean(common$se_unadj^2)),
            CP = mean(s$covered), mean_CI_width = mean(s$ci_high -
                s$ci_low), reject_rate = mean(s$reject), bias_mcse = sd(s$estimate)/sqrt(nrow(s)),
            max_gram_condition = if (nrow(s)) max(s$gram_condition) else NA_real_,
            max_normal_error = if (nrow(s)) max(s$normal_error) else NA_real_)
    }, by = .(N, effect, W, estimand)]
    out[, `:=`(type_I_error = ifelse(effect == "null", reject_rate,
        NA_real_), power = ifelse(effect == "alternative", reject_rate,
        NA_real_))]
    ci <- wilson_ci(out$reject_rate, out$n_success)
    out[, `:=`(reject_mc_low = ci[, 1], reject_mc_high = ci[,
        2])]
    ci <- wilson_ci(out$CP, out$n_success)
    out[, `:=`(CP_mc_low = ci[, 1], CP_mc_high = ci[, 2])]
    out <- merge(out, oracle[, .(effect, W, estimand, oracle_PVR,
        conditional_mcse_PVR)], by = c("effect", "W", "estimand"),
        all.x = TRUE)
    out[, `:=`(PVR_minus_oracle, empirical_PVR - oracle_PVR)]
    data.table::setorder(out, effect, N, estimand, W)
    save_csv(out, "main_simulation_results", output_dir)
    failures <- raw[success == FALSE, .N, by = .(N, effect, W,
        estimand, failure_reason)]
    save_csv(failures, "numerical_failures", output_dir)
    save_csv(raw[success == FALSE], "failed_replications", output_dir)
    flags <- out[, .(N, effect, W, estimand, n_success, n_failed,
        ASE_ESD, RMSSE_ESD, CP, CP_mc_low, CP_mc_high, reject_rate,
        reject_mc_low, reject_mc_high, empirical_PVR, oracle_PVR,
        PVR_minus_oracle, SE_discrepancy = abs(RMSSE_ESD - 1) >
            0.1, coverage_nominal_outside_MC_CI = CP_mc_high <
            0.95 | CP_mc_low > 0.95, type1_nominal_outside_MC_CI = effect ==
            "null" & (reject_mc_high < 0.05 | reject_mc_low >
            0.05), finite_sample_variance_loss = empirical_PVR <
            0)]
    save_csv(flags, "operating_characteristic_diagnostics", output_dir)
    out
}
