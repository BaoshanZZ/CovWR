local_project_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 1L) {
    return(dirname(dirname(normalizePath(sub("^--file=", "", file_arg)))))
  }
  if (dir.exists("WRCov/R")) {
    return(normalizePath("WRCov"))
  }
  normalizePath(".")
}

parse_validation_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(hit) == 1L) {
    sub(paste0("^--", name, "="), "", hit)
  } else {
    default
  }
}

project_dir <- local_project_dir()
r_dir <- file.path(project_dir, "R")
source(file.path(r_dir, "config.R"))
source(file.path(r_dir, "dgm.R"))
source(file.path(r_dir, "win_kernels.R"))
source(file.path(r_dir, "estimators.R"))

mc_reps <- as.integer(parse_validation_arg("mc-reps", "500"))
n_cores <- as.integer(parse_validation_arg("cores", "1"))
output_dir <- parse_validation_arg(
  "output-dir",
  file.path(project_dir, "formula_validation")
)
if (is.na(mc_reps) || mc_reps < 1L) {
  stop("--mc-reps must be a positive integer.", call. = FALSE)
}
if (is.na(n_cores) || n_cores < 1L) {
  stop("--cores must be a positive integer.", call. = FALSE)
}

checks <- list()
add_check <- function(name, observed, target, tolerance, relation = "absolute") {
  error <- switch(
    relation,
    absolute = abs(observed - target),
    less_equal = observed - target,
    greater_equal = target - observed,
    stop("Unknown relation.", call. = FALSE)
  )
  passed <- is.finite(error) && error <= tolerance
  checks[[length(checks) + 1L]] <<- data.table::data.table(
    check = name,
    observed = observed,
    target = target,
    tolerance = tolerance,
    relation = relation,
    error = error,
    passed = passed
  )
  invisible(passed)
}

params <- default_dgm_params()
set.seed(2026072601)
dat <- generate_trial_data(
  N = 240L,
  e = 0.4,
  effect = "null",
  prognosis = "moderate",
  params = params
)
pairwise <- compute_pairwise_quantities(
  dat = dat,
  allocation_prob = 0.4,
  params = params
)

add_check(
  "Hierarchical kernel partitions every pair into win/loss/tie",
  max(abs(pairwise$matrices$win + pairwise$matrices$loss + pairwise$matrices$tie - 1)),
  0,
  0
)
add_check(
  "Win and loss kernels are mutually exclusive",
  max(pairwise$matrices$win * pairwise$matrices$loss),
  0,
  0
)

# Exact Hoeffding-Hajek decomposition for a finite-support win kernel.
treated_binary <- c(0, 1, 1, 0, 1)
control_binary <- c(1, 0, 0, 1)
binary_kernel <- outer(treated_binary, control_binary, ">") * 1
binary_tau <- 0.25
binary_g_t <- ifelse(treated_binary == 1, 0.25, -0.25)
binary_g_c <- ifelse(control_binary == 0, 0.25, -0.25)
binary_projection <- mean(binary_g_t) + mean(binary_g_c)
binary_residual <- binary_kernel -
  binary_tau -
  outer(binary_g_t, binary_g_c, "+")
add_check(
  "Hoeffding-Hajek projection and degenerate remainder decompose U exactly",
  mean(binary_kernel) - binary_tau,
  binary_projection + mean(binary_residual),
  1e-15
)

boundary_treated <- data.table::data.table(
  death_time = params$S,
  death_event = 0L,
  alsfrs = 38,
  nfl_score = params$nfl_epsilon
)
boundary_control <- data.table::data.table(
  death_time = params$S,
  death_event = 0L,
  alsfrs = 36,
  nfl_score = 0
)
boundary_pair <- hierarchical_pair_matrices(boundary_treated, boundary_control, params)
add_check(
  "Differences exactly equal to continuous thresholds remain tied",
  boundary_pair$tie[1L, 1L],
  1,
  0
)

add_check(
  "U-statistic win probability equals the pairwise mean",
  pairwise$tau_w,
  mean(pairwise$matrices$win),
  1e-15
)
add_check(
  "U-statistic loss probability equals the pairwise mean",
  pairwise$tau_l,
  mean(pairwise$matrices$loss),
  1e-15
)

treated_idx <- dat$A == 1L
control_idx <- !treated_idx
placement <- pairwise$placement
placement_centering_error <- max(abs(c(
  mean(placement$g_w_T[treated_idx]),
  mean(placement$g_l_T[treated_idx]),
  mean(placement$g_w_C[control_idx]),
  mean(placement$g_l_C[control_idx])
)))
add_check(
  "All four empirical placement functions are arm-wise centered",
  placement_centering_error,
  0,
  1e-14
)
add_check(
  "Estimated elementary IFs have empirical mean zero",
  max(abs(colMeans(pairwise$psi))),
  0,
  1e-14
)
add_check(
  "Elementary IF uses the prespecified allocation probability",
  pairwise$allocation_prob,
  0.4,
  0
)

tau_for_gradient <- c(w = 0.43, l = 0.31)
step <- 1e-6
for (measure in c("NB", "LWR", "LWO")) {
  analytic <- measure_gradient(
    measure,
    tau_for_gradient[["w"]],
    tau_for_gradient[["l"]],
    params
  )
  numeric_gradient <- c(
    w = (
      calc_measures(tau_for_gradient[["w"]] + step, tau_for_gradient[["l"]], params)[[measure]] -
        calc_measures(tau_for_gradient[["w"]] - step, tau_for_gradient[["l"]], params)[[measure]]
    ) / (2 * step),
    l = (
      calc_measures(tau_for_gradient[["w"]], tau_for_gradient[["l"]] + step, params)[[measure]] -
        calc_measures(tau_for_gradient[["w"]], tau_for_gradient[["l"]] - step, params)[[measure]]
    ) / (2 * step)
  )
  add_check(
    paste(measure, "analytic gradient matches numerical differentiation"),
    max(abs(analytic - numeric_gradient)),
    0,
    1e-8
  )
}

N <- nrow(dat)
for (measure in c("NB", "LWR", "LWO")) {
  grad <- measure_gradient(measure, pairwise$tau_w, pairwise$tau_l, params)
  score <- as.vector(pairwise$psi %*% grad)
  direct_var <- sum((score - mean(score))^2) / N^2
  covariance_var <- if_se_for_measure(
    pairwise$psi,
    pairwise$tau_w,
    pairwise$tau_l,
    measure,
    N,
    params
  )^2
  add_check(
    paste(measure, "unadjusted covariance and score variance formulas agree"),
    covariance_var,
    direct_var,
    1e-15
  )
}

ladj <- fit_if_projection_components(dat, pairwise)
Fmat <- adjustment_features(dat)
Z <- dat$A - pairwise$allocation_prob
design <- Fmat * as.numeric(Z)
add_check(
  "LADJ feature columns are sample centered",
  max(abs(colMeans(Fmat))),
  0,
  1e-14
)
add_check(
  "LADJ uses no intercept column",
  as.numeric("intercept" %in% colnames(Fmat)),
  0,
  0
)
normal_equation_error <- max(abs(c(
  crossprod(design, pairwise$psi[, "w"] - design %*% ladj$beta[, "w"]),
  crossprod(design, pairwise$psi[, "l"] - design %*% ladj$beta[, "l"])
)))
add_check(
  "Closed-form LADJ coefficients satisfy both normal equations",
  normal_equation_error,
  0,
  1e-10
)

mean_design <- colMeans(design)
centering_identity <- ladj$imbalance_factor * ladj$imbalance
add_check(
  "Sample-centering imbalance identity holds",
  max(abs(mean_design - centering_identity)),
  0,
  1e-14
)
tau_direct <- c(
  w = pairwise$tau_w - mean(Z * ladj$b_hat[, "w"]),
  l = pairwise$tau_l - mean(Z * ladj$b_hat[, "l"])
)
add_check(
  "LADJ imbalance and direct augmentation estimators agree",
  max(abs(tau_direct - ladj$tau_adj)),
  0,
  1e-14
)
direct_risk <- sum((pairwise$psi[, "w"] - Z * ladj$b_hat[, "w"])^2)
weighted_risk <- sum(Z^2 * (pairwise$psi[, "w"] / Z - ladj$b_hat[, "w"])^2)
add_check(
  "Direct and weighted flexible-learning risk criteria are identical",
  direct_risk,
  weighted_risk,
  1e-12
)
add_check(
  "Adjusted elementary IF equals psi minus (A-e)b(X)",
  max(abs(ladj$psi_adj - (pairwise$psi - Z * ladj$b_hat))),
  0,
  0
)

for (measure in c("NB", "LWR", "LWO")) {
  grad <- measure_gradient(measure, pairwise$tau_w, pairwise$tau_l, params)
  beta_combined <- as.vector(ladj$beta %*% grad)
  beta_direct <- as.vector(solve(crossprod(design), crossprod(
    design,
    as.vector(pairwise$psi %*% grad)
  )))
  add_check(
    paste(measure, "common elementary projection equals direct measure projection"),
    max(abs(beta_combined - beta_direct)),
    0,
    1e-12
  )

  unadjusted_score <- as.vector(pairwise$psi %*% grad)
  adjusted_score <- as.vector(ladj$psi_adj %*% grad)
  empirical_unadjusted_variance <- mean(
    (unadjusted_score - mean(unadjusted_score))^2
  )
  empirical_adjusted_variance <- mean(
    (adjusted_score - mean(adjusted_score))^2
  )
  add_check(
    paste(measure, "linear projection has finite-sample empirical no-loss"),
    empirical_adjusted_variance,
    empirical_unadjusted_variance,
    1e-12,
    relation = "less_equal"
  )

  adjusted_covariance_var <- if_se_for_measure(
    ladj$psi_adj,
    ladj$tau_adj[["w"]],
    ladj$tau_adj[["l"]],
    measure,
    N,
    params
  )^2
  adjusted_grad <- measure_gradient(
    measure,
    ladj$tau_adj[["w"]],
    ladj$tau_adj[["l"]],
    params
  )
  adjusted_measure_score <- as.vector(ladj$psi_adj %*% adjusted_grad)
  adjusted_score_var <- sum(
    (adjusted_measure_score - mean(adjusted_measure_score))^2
  ) / N^2
  add_check(
    paste(measure, "adjusted covariance and score variance formulas agree"),
    adjusted_covariance_var,
    adjusted_score_var,
    1e-15
  )
}

null_gradients <- rbind(
  NB = measure_gradient("NB", 0.4, 0.4, params),
  LWR = measure_gradient("LWR", 0.4, 0.4, params),
  LWO = measure_gradient("LWO", 0.4, 0.4, params)
)
gradient_ratios <- null_gradients[, "w"] / c(NB = 1, LWR = 2.5, LWO = 2)
add_check(
  "All three null gradients are proportional to (1,-1)",
  max(abs(gradient_ratios - 1), abs(null_gradients[, "w"] + null_gradients[, "l"])),
  0,
  1e-14
)

# Exact finite-support check of the optimal augmentation formula.
e_opt <- 0.4
support <- data.table::CJ(X = c(-1, 1), A = 0:1, R = c(-1, 1))
support[, probability := 0.5 *
  ifelse(A == 1L, e_opt, 1 - e_opt) *
  0.5]
support[, psi := ifelse(A == 1L, 0.3 + 0.8 * X, -0.2 - 0.4 * X) + 0.2 * R]
support[, b_opt := (0.3 + 0.8 * X) - (-0.2 - 0.4 * X)]
support[, augmentation := (A - e_opt) * b_opt]
support[, residual := psi - augmentation]
weighted_mean <- function(x) sum(support$probability * x)
weighted_var <- function(x) {
  center <- weighted_mean(x)
  weighted_mean((x - center)^2)
}
orthogonality <- weighted_mean(support$residual * support$augmentation)
pythagorean_error <- weighted_var(support$psi) -
  weighted_var(support$residual) -
  weighted_var(support$augmentation)
add_check(
  "Optimal conditional-contrast residual is orthogonal to augmentation",
  orthogonality,
  0,
  1e-14
)
add_check(
  "Randomization makes the optimal augmentation term mean zero",
  weighted_mean(support$augmentation),
  0,
  1e-14
)
add_check(
  "Optimal augmentation Pythagorean variance identity holds",
  pythagorean_error,
  0,
  1e-14
)
add_check(
  "Optimal augmentation has no larger risk than zero augmentation",
  weighted_mean(support$residual^2),
  weighted_mean(support$psi^2),
  1e-14,
  relation = "less_equal"
)

# Exact finite-support check of the centered linear variance reduction formula.
linear_support <- data.table::CJ(X = c(-1, 1), A = 0:1, R = c(-1, 1))
linear_support[, probability := 0.5 *
  ifelse(A == 1L, e_opt, 1 - e_opt) *
  0.5]
linear_support[, psi := ifelse(A == 1L, 0.8 * X, -0.4 * X) + 0.2 * R]
linear_support[, design := (A - e_opt) * X]
linear_weighted_mean <- function(x) sum(linear_support$probability * x)
linear_weighted_var <- function(x) {
  center <- linear_weighted_mean(x)
  linear_weighted_mean((x - center)^2)
}
gamma_linear <- linear_weighted_mean(linear_support$design * linear_support$psi) /
  linear_weighted_mean(linear_support$design^2)
linear_support[, residual := psi - design * gamma_linear]
linear_variance_reduction <- linear_weighted_var(linear_support$psi) -
  linear_weighted_var(linear_support$residual)
linear_variance_formula <- e_opt * (1 - e_opt) *
  gamma_linear^2 *
  linear_weighted_var(linear_support$X)
add_check(
  "Arm-wise centered IF has zero population intercept coefficient",
  linear_weighted_mean((linear_support$A - e_opt) * linear_support$psi),
  0,
  1e-14
)
add_check(
  "Centered linear projection recovers the conditional-contrast slope",
  gamma_linear,
  1.2,
  1e-14
)
add_check(
  "Centered linear variance-reduction closed form holds",
  linear_variance_reduction,
  linear_variance_formula,
  1e-14
)

exact_results <- data.table::rbindlist(checks)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(
  exact_results,
  file.path(output_dir, "exact_formula_checks.csv")
)

run_mc_replication <- function(N, rep_id) {
  set.seed(2026072700 + N * 10000L + rep_id)
  trial <- generate_trial_data(
    N = N,
    e = 0.5,
    effect = "null",
    prognosis = "moderate",
    params = params
  )
  estimate_all_methods(
    dat = trial,
    allocation_prob = 0.5,
    truth_values = c(NB = 0, LWR = 0, LWO = 0),
    params = params,
    method_names = c("UNADJ", "LADJ")
  )[, .(
    N = N,
    rep = rep_id,
    method,
    estimand,
    estimate,
    se,
    reject
  )]
}

mc_grid <- data.table::CJ(N = c(300L, 1000L), rep = seq_len(mc_reps))
run_mc_task <- function(k) {
  run_mc_replication(mc_grid$N[[k]], mc_grid$rep[[k]])
}
if (n_cores > 1L && .Platform$OS.type != "windows") {
  mc_rows <- parallel::mclapply(
    seq_len(nrow(mc_grid)),
    run_mc_task,
    mc.cores = n_cores,
    mc.preschedule = FALSE
  )
} else {
  mc_rows <- lapply(seq_len(nrow(mc_grid)), run_mc_task)
}
mc_raw <- data.table::rbindlist(mc_rows)
mc_summary <- mc_raw[, {
  esd <- stats::sd(estimate)
  rmsse <- sqrt(mean(se^2))
  z <- estimate / se
  reject_rate <- mean(reject)
  mc_se_reject <- sqrt(reject_rate * (1 - reject_rate) / .N)
  .(
    MC_reps = .N,
    mean_estimate = mean(estimate),
    ESD = esd,
    RMSSE = rmsse,
    RMSSE_ESD = rmsse / esd,
    mean_z = mean(z),
    SD_z = stats::sd(z),
    reject_rate = reject_rate,
    reject_mc_low = max(0, reject_rate - stats::qnorm(0.975) * mc_se_reject),
    reject_mc_high = min(1, reject_rate + stats::qnorm(0.975) * mc_se_reject)
  )
}, by = .(N, method, estimand)]
data.table::fwrite(
  mc_summary,
  file.path(output_dir, "asymptotic_mc_checks.csv")
)

cat("\nExact formula checks:\n")
print(exact_results[, .(check, passed, error, tolerance)])
cat("\nMonte Carlo asymptotic checks:\n")
print(mc_summary)
cat("\nSaved validation output to:", normalizePath(output_dir), "\n")

if (!all(exact_results$passed)) {
  stop("At least one exact formula check failed.", call. = FALSE)
}
