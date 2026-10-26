# Mechanically extracted from the validated retained implementation; estimator logic is unchanged.
legacy_survivor_data <- function (dat)
{
    data.table::data.table(death_time = dat$death_time, death_event = dat$death_event,
        alsfrs = -dat$icu_days, nfl_score = -dat$log_biomarker)
}

compare_hce_pair <- function (treated, control, params)
{
    stopifnot(nrow(treated) == nrow(control))
    win <- treated$death_time > control$death_time & control$death_event ==
        1
    loss <- treated$death_time < control$death_time & treated$death_event ==
        1
    tier <- ifelse(win | loss, "D1", "overall tie")
    survivors <- which(treated$death_event == 0 & control$death_event ==
        0)
    if (length(survivors)) {
        t <- treated[survivors]
        c <- control[survivors]
        stopifnot(!anyNA(t$log_biomarker), !anyNA(c$log_biomarker))
        legacy <- legacy_pair_vectors(legacy_survivor_data(t),
            legacy_survivor_data(c), legacy_margins)
        win[survivors] <- as.logical(legacy$win)
        loss[survivors] <- as.logical(legacy$loss)
        resolved2 <- abs(t$icu_days - c$icu_days) > params$icu_margin
        tier[survivors] <- ifelse(resolved2, "D2", ifelse(legacy$win |
            legacy$loss, "D3", "overall tie"))
    }
    data.table::data.table(win = as.integer(win), loss = as.integer(loss),
        tie = as.integer(!(win | loss)), deciding_component = tier)
}

revised_pair_matrices <- function (treated, control, params)
{
    n1 <- nrow(treated)
    n0 <- nrow(control)
    if (!n1 || !n0)
        stop("One randomization arm is empty.")
    win <- outer(treated$death_time, control$death_time, ">") &
        matrix(control$death_event == 1, n1, n0, byrow = TRUE)
    loss <- outer(treated$death_time, control$death_time, "<") &
        matrix(treated$death_event == 1, n1, n0)
    tier <- matrix(ifelse(win | loss, 1L, 0L), n1, n0)
    ti <- which(treated$death_event == 0)
    ci <- which(control$death_event == 0)
    if (length(ti) && length(ci)) {
        t <- treated[ti]
        c <- control[ci]
        stopifnot(!anyNA(t$log_biomarker), !anyNA(c$log_biomarker))
        legacy <- legacy_pair_matrices(legacy_survivor_data(t),
            legacy_survivor_data(c), legacy_margins)
        win[ti, ci] <- legacy$win
        loss[ti, ci] <- legacy$loss
        d2 <- abs(outer(t$icu_days, c$icu_days, "-")) > params$icu_margin
        tier[ti, ci] <- ifelse(d2, 2L, ifelse(legacy$win | legacy$loss,
            3L, 0L))
    }
    list(win = win * 1, loss = loss * 1, tie = (!(win | loss)) *
        1, deciding_component = tier)
}

make_reference_bank <- function (reference, params)
{
    alive <- reference[death_event == 0]
    list(n = nrow(reference), times = sort(reference$death_time),
        deaths = sort(reference$death_time[reference$death_event ==
            1]), icu = sort(alive$icu_days), scores = lapply(0:params$S,
            function(k) sort(-alive$log_biomarker[alive$icu_days ==
                k])))
}

bank_pair_expectations <- function (target, bank, params)
{
    less <- function(x, sorted) findInterval(x, sorted, left.open = TRUE)
    leq <- function(x, sorted) findInterval(x, sorted)
    win <- less(target$death_time, bank$deaths)
    loss <- ifelse(target$death_event == 1, bank$n - leq(target$death_time,
        bank$times), 0)
    ti <- which(target$death_event == 0)
    if (length(ti)) {
        x <- target$icu_days[ti]
        score <- -target$log_biomarker[ti]
        stopifnot(!anyNA(score))
        win[ti] <- win[ti] + length(bank$icu) - leq(x + params$icu_margin,
            bank$icu)
        loss[ti] <- loss[ti] + less(x - params$icu_margin, bank$icu)
        for (d in 0:params$S) {
            use <- which(abs(x - d) <= params$icu_margin)
            ref <- bank$scores[[d + 1L]]
            if (!length(use) || !length(ref))
                next
            win[ti[use]] <- win[ti[use]] + less(score[use], ref +
                params$biomarker_margin)
            loss[ti[use]] <- loss[ti[use]] + length(ref) - leq(score[use],
                ref - params$biomarker_margin)
        }
    }
    cbind(w = win/bank$n, l = loss/bank$n)
}

population_placements <- function (target, reference, params, block_size = 10000L)
{
    bank <- make_reference_bank(reference, params)
    result <- matrix(NA_real_, nrow(target), 2, dimnames = list(NULL,
        c("w", "l")))
    for (start in seq.int(1L, nrow(target), by = block_size)) {
        idx <- start:min(nrow(target), start + block_size - 1L)
        result[idx, ] <- bank_pair_expectations(target[idx],
            bank, params)
    }
    result
}
