hierarchical_pair_matrices <- function(treated, control, params = default_dgm_params()) {
  n1 <- nrow(treated)
  n0 <- nrow(control)
  if (n1 == 0L || n0 == 0L) {
    stop("Both treatment arms must contain at least one subject.", call. = FALSE)
  }

  death_win <- outer(treated$death_time, control$death_time + params$death_epsilon, ">") &
    matrix(control$death_event == 1L, nrow = n1, ncol = n0, byrow = TRUE)
  death_loss <- outer(treated$death_time, control$death_time - params$death_epsilon, "<") &
    matrix(treated$death_event == 1L, nrow = n1, ncol = n0)
  death_tie <- !(death_win | death_loss)

  als_win <- outer(treated$alsfrs, control$alsfrs + params$alsfrs_epsilon, ">")
  als_loss <- outer(treated$alsfrs, control$alsfrs - params$alsfrs_epsilon, "<")
  als_tie <- !(als_win | als_loss)

  nfl_win <- outer(treated$nfl_score, control$nfl_score + params$nfl_epsilon, ">")
  nfl_loss <- outer(treated$nfl_score, control$nfl_score - params$nfl_epsilon, "<")

  win <- death_win | (death_tie & als_win) | (death_tie & als_tie & nfl_win)
  loss <- death_loss | (death_tie & als_loss) | (death_tie & als_tie & nfl_loss)
  tie <- !(win | loss)

  list(
    win = matrix(as.numeric(win), nrow = n1, ncol = n0),
    loss = matrix(as.numeric(loss), nrow = n1, ncol = n0),
    tie = matrix(as.numeric(tie), nrow = n1, ncol = n0)
  )
}

hierarchical_pair_vectors <- function(treated, control, params = default_dgm_params()) {
  if (nrow(treated) != nrow(control)) {
    stop("treated and control must have the same number of rows for paired truth estimation.", call. = FALSE)
  }

  death_win <- (treated$death_time > control$death_time + params$death_epsilon) &
    (control$death_event == 1L)
  death_loss <- (treated$death_time < control$death_time - params$death_epsilon) &
    (treated$death_event == 1L)
  death_tie <- !(death_win | death_loss)

  als_win <- treated$alsfrs > control$alsfrs + params$alsfrs_epsilon
  als_loss <- treated$alsfrs < control$alsfrs - params$alsfrs_epsilon
  als_tie <- !(als_win | als_loss)

  nfl_win <- treated$nfl_score > control$nfl_score + params$nfl_epsilon
  nfl_loss <- treated$nfl_score < control$nfl_score - params$nfl_epsilon

  win <- death_win | (death_tie & als_win) | (death_tie & als_tie & nfl_win)
  loss <- death_loss | (death_tie & als_loss) | (death_tie & als_tie & nfl_loss)
  tie <- !(win | loss)

  data.table::data.table(
    win = as.numeric(win),
    loss = as.numeric(loss),
    tie = as.numeric(tie)
  )
}

check_kernel_partition <- function(pair_obj) {
  sums <- pair_obj$win + pair_obj$loss + pair_obj$tie
  isTRUE(all.equal(as.numeric(sums), rep(1, length(sums)), tolerance = 0))
}

compute_pairwise_quantities <- function(dat, params = default_dgm_params()) {
  treated <- dat[A == 1L]
  control <- dat[A == 0L]
  mats <- hierarchical_pair_matrices(treated, control, params)

  tau_w <- mean(mats$win)
  tau_l <- mean(mats$loss)

  g_w_T <- rowMeans(mats$win) - tau_w
  g_w_C <- colMeans(mats$win) - tau_w
  g_l_T <- rowMeans(mats$loss) - tau_l
  g_l_C <- colMeans(mats$loss) - tau_l

  placement <- data.table::data.table(
    id = c(treated$id, control$id),
    g_w_T = c(g_w_T, rep(NA_real_, nrow(control))),
    g_w_C = c(rep(NA_real_, nrow(treated)), g_w_C),
    g_l_T = c(g_l_T, rep(NA_real_, nrow(control))),
    g_l_C = c(rep(NA_real_, nrow(treated)), g_l_C)
  )
  data.table::setkey(placement, id)
  placement <- placement[dat[, .(id)]]

  ehat <- mean(dat$A)
  zero_na <- function(x) {
    x[is.na(x)] <- 0
    x
  }
  psi_w <- dat$A / ehat * zero_na(placement$g_w_T) +
    (1 - dat$A) / (1 - ehat) * zero_na(placement$g_w_C)
  psi_l <- dat$A / ehat * zero_na(placement$g_l_T) +
    (1 - dat$A) / (1 - ehat) * zero_na(placement$g_l_C)

  list(
    tau_w = tau_w,
    tau_l = tau_l,
    ehat = ehat,
    placement = placement,
    psi = cbind(w = psi_w, l = psi_l),
    matrices = mats
  )
}
