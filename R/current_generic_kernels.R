# Mechanically extracted from the validated retained implementation; estimator logic is unchanged.
legacy_pair_matrices <- function (treated, control, params = default_dgm_params())
{
    n1 <- nrow(treated)
    n0 <- nrow(control)
    if (n1 == 0L || n0 == 0L) {
        stop("Both treatment arms must contain at least one subject.",
            call. = FALSE)
    }
    death_win <- outer(treated$death_time, control$death_time +
        params$death_epsilon, ">") & matrix(control$death_event ==
        1L, nrow = n1, ncol = n0, byrow = TRUE)
    death_loss <- outer(treated$death_time, control$death_time -
        params$death_epsilon, "<") & matrix(treated$death_event ==
        1L, nrow = n1, ncol = n0)
    death_tie <- !(death_win | death_loss)
    als_win <- outer(treated$alsfrs, control$alsfrs + params$alsfrs_epsilon,
        ">")
    als_loss <- outer(treated$alsfrs, control$alsfrs - params$alsfrs_epsilon,
        "<")
    als_tie <- !(als_win | als_loss)
    nfl_win <- outer(treated$nfl_score, control$nfl_score + params$nfl_epsilon,
        ">")
    nfl_loss <- outer(treated$nfl_score, control$nfl_score -
        params$nfl_epsilon, "<")
    win <- death_win | (death_tie & als_win) | (death_tie & als_tie &
        nfl_win)
    loss <- death_loss | (death_tie & als_loss) | (death_tie &
        als_tie & nfl_loss)
    tie <- !(win | loss)
    list(win = matrix(as.numeric(win), nrow = n1, ncol = n0),
        loss = matrix(as.numeric(loss), nrow = n1, ncol = n0),
        tie = matrix(as.numeric(tie), nrow = n1, ncol = n0))
}

legacy_pair_vectors <- function (treated, control, params = default_dgm_params())
{
    if (nrow(treated) != nrow(control)) {
        stop("treated and control must have the same number of rows for paired truth estimation.",
            call. = FALSE)
    }
    death_win <- (treated$death_time > control$death_time + params$death_epsilon) &
        (control$death_event == 1L)
    death_loss <- (treated$death_time < control$death_time -
        params$death_epsilon) & (treated$death_event == 1L)
    death_tie <- !(death_win | death_loss)
    als_win <- treated$alsfrs > control$alsfrs + params$alsfrs_epsilon
    als_loss <- treated$alsfrs < control$alsfrs - params$alsfrs_epsilon
    als_tie <- !(als_win | als_loss)
    nfl_win <- treated$nfl_score > control$nfl_score + params$nfl_epsilon
    nfl_loss <- treated$nfl_score < control$nfl_score - params$nfl_epsilon
    win <- death_win | (death_tie & als_win) | (death_tie & als_tie &
        nfl_win)
    loss <- death_loss | (death_tie & als_loss) | (death_tie &
        als_tie & nfl_loss)
    tie <- !(win | loss)
    data.table::data.table(win = as.numeric(win), loss = as.numeric(loss),
        tie = as.numeric(tie))
}
