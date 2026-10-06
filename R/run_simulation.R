local_r_dir <- function() {
  frames <- sys.frames()
  files <- vapply(frames, function(frame) {
    if (!is.null(frame$ofile)) {
      frame$ofile
    } else {
      NA_character_
    }
  }, character(1))
  files <- stats::na.omit(files)
  if (length(files) > 0L) {
    file_dir <- dirname(normalizePath(files[[length(files)]]))
    if (file.exists(file.path(file_dir, "config.R"))) {
      return(file_dir)
    }
  }
  if (dir.exists("WRCov/R")) {
    return(normalizePath("WRCov/R"))
  }
  if (dir.exists("R") && file.exists(file.path("R", "config.R"))) {
    return(normalizePath("R"))
  }
  normalizePath(".")
}

.wrcov_r_dir <- local_r_dir()
source(file.path(.wrcov_r_dir, "config.R"))
check_sim_packages()
source(file.path(.wrcov_r_dir, "dgm.R"))
source(file.path(.wrcov_r_dir, "win_kernels.R"))
source(file.path(.wrcov_r_dir, "estimators.R"))
source(file.path(.wrcov_r_dir, "truth.R"))
source(file.path(.wrcov_r_dir, "summarize_results.R"))

run_one_replication <- function(scenario,
                                rep_id,
                                rep_seed,
                                truth_row,
                                params,
                                sbw_nonnegative,
                                bootstrap_B = 0L,
                                method_names = c("UNADJ", "LADJ")) {
  set.seed(rep_seed)
  dat <- generate_trial_data(
    N = scenario$N,
    e = scenario$e,
    effect = scenario$effect,
    prognosis = scenario$prognosis,
    outcome_shape = scenario$outcome_shape,
    effect_heterogeneity = scenario$effect_heterogeneity,
    endpoint_prognosis_pattern = scenario$endpoint_prognosis_pattern,
    n_noise_covariates = scenario$n_noise_covariates,
    adjustment_covariate_pattern = scenario$adjustment_covariate_pattern,
    params = params
  )

  truth_values <- truth_vector_for_scenario(truth_row)
  est <- estimate_all_methods(
    dat = dat,
    allocation_prob = scenario$e,
    truth_values = truth_values,
    params = params,
    sbw_nonnegative = sbw_nonnegative,
    bootstrap_B = bootstrap_B,
    method_names = method_names
  )

  est[, `:=`(
    scenario_id = scenario$scenario_id,
    rep = rep_id,
    N = scenario$N,
    e = scenario$e,
    effect = scenario$effect,
    prognosis = scenario$prognosis,
    outcome_shape = scenario$outcome_shape,
    effect_heterogeneity = scenario$effect_heterogeneity,
    endpoint_prognosis_pattern = scenario$endpoint_prognosis_pattern,
    n_noise_covariates = scenario$n_noise_covariates,
    adjustment_covariate_pattern = scenario$adjustment_covariate_pattern
  )]

  data.table::setcolorder(est, c(
    "scenario_id", "rep", "N", "e", "effect", "prognosis",
    "outcome_shape", "effect_heterogeneity", "endpoint_prognosis_pattern",
    "n_noise_covariates", "adjustment_covariate_pattern",
    "method", "estimand", "estimate", "se", "ci_low", "ci_high", "p_value",
    "reject", "truth", "covered", "tau_w", "tau_l", "learner_w", "learner_l",
    "bootstrap_B", "boot_var", "se_boot", "ci_low_boot", "ci_high_boot",
    "p_value_boot", "reject_boot", "covered_boot"
  ))
  est[]
}

run_simulation <- function(run_full = FALSE,
                           n_cores = 1L,
                           seed = 20260617,
                           scenario_grid = NULL,
                           scenario_design = c("one_factor", "full_factorial"),
                           bootstrap_B = NULL,
                           method_names = c("UNADJ", "LADJ"),
                           n_rep = NULL,
                           rep_start = 0L,
                           truth_override = NULL,
                           write_outputs = FALSE,
                           output_dir = file.path(dirname(.wrcov_r_dir), "simulation_output")) {
  params <- default_dgm_params()
  settings <- simulation_settings(run_full)
  if (!is.null(bootstrap_B)) {
    settings$bootstrap_B <- as.integer(bootstrap_B)
  }
  if (!is.null(n_rep)) {
    settings$n_rep <- as.integer(n_rep)
  }
  rep_start <- as.integer(rep_start)
  scenario_design <- match.arg(scenario_design)
  grid <- if (is.null(scenario_grid)) {
    make_scenario_grid(run_full, scenario_design = scenario_design)
  } else {
    normalize_scenario_grid(scenario_grid)
  }

  if (is.null(truth_override)) {
    message("Estimating truth values for ", nrow(grid), " scenario(s).")
    truth <- estimate_truth_grid(
      grid = grid,
      n_pairs = settings$truth_pairs,
      chunk_size = settings$truth_chunk_size,
      seed = seed,
      params = params
    )
  } else {
    truth <- data.table::as.data.table(truth_override)
    required_truth_cols <- c(
      "scenario_id", "truth_tau_w", "truth_tau_l", "truth_NB", "truth_LWR", "truth_LWO"
    )
    missing_truth_cols <- setdiff(required_truth_cols, names(truth))
    if (length(missing_truth_cols) > 0L) {
      stop("truth_override is missing column(s): ", paste(missing_truth_cols, collapse = ", "), call. = FALSE)
    }
    truth <- truth[scenario_id %in% grid$scenario_id, ..required_truth_cols]
    if (!setequal(truth$scenario_id, grid$scenario_id)) {
      stop("truth_override must contain exactly one row for each scenario_id.", call. = FALSE)
    }
    message("Using supplied truth values for ", nrow(grid), " scenario(s).")
  }

  tasks <- data.table::CJ(
    scenario_id = grid$scenario_id,
    rep = seq.int(rep_start + 1L, rep_start + settings$n_rep),
    sorted = FALSE
  )
  rng_id <- stats::setNames(grid$rng_scenario_id, grid$scenario_id)
  tasks[, rep_seed := seed + rng_id[as.character(scenario_id)] * 100000L + rep]

  run_task <- function(k) {
    task <- tasks[k, , drop = FALSE]
    sid <- task[["scenario_id"]]
    scenario <- grid[grid[["scenario_id"]] == sid, , drop = FALSE]
    truth_row <- truth[truth[["scenario_id"]] == sid, , drop = FALSE]
    run_one_replication(
      scenario = scenario,
      rep_id = task[["rep"]],
      rep_seed = task[["rep_seed"]],
      truth_row = truth_row,
      params = params,
      sbw_nonnegative = settings$sbw_nonnegative,
      bootstrap_B = settings$bootstrap_B,
      method_names = method_names
    )
  }

  message("Running ", nrow(tasks), " replication task(s).")
  if (n_cores > 1L) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n_cores)

    chunk_size <- max(1L, min(100L, ceiling(nrow(tasks) / (n_cores * 20L))))
    task_chunks <- split(seq_len(nrow(tasks)), ceiling(seq_len(nrow(tasks)) / chunk_size))
    run_chunk <- function(idx) {
      data.table::rbindlist(lapply(idx, run_task), use.names = TRUE)
    }
    raw <- data.table::rbindlist(
      future.apply::future_lapply(task_chunks, run_chunk, future.seed = TRUE),
      use.names = TRUE
    )
  } else {
    raw <- data.table::rbindlist(lapply(seq_len(nrow(tasks)), run_task), use.names = TRUE)
  }

  summary <- summarize_simulation_results(raw)
  check_simulation_outputs(raw, summary)

  result <- list(
    raw = raw,
    summary = summary,
    truth = truth,
    scenario_grid = grid,
    settings = settings,
    params = params
  )

  if (write_outputs) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(raw, file.path(output_dir, "simulation_raw.csv"))
    data.table::fwrite(summary, file.path(output_dir, "simulation_summary.csv"))
    data.table::fwrite(truth, file.path(output_dir, "simulation_truth.csv"))
  }

  result
}

running_as_direct_script <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) != 1L) {
    return(FALSE)
  }
  script_file <- normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE)
  this_file <- normalizePath(file.path(.wrcov_r_dir, "run_simulation.R"), mustWork = FALSE)
  identical(script_file, this_file)
}

if (running_as_direct_script()) {
  args <- commandArgs(trailingOnly = TRUE)
  run_full <- "--full" %in% args
  write_outputs <- "--write" %in% args
  n_cores_arg <- grep("^--cores=", args, value = TRUE)
  bootstrap_arg <- grep("^--bootstrap-B=", args, value = TRUE)
  methods_arg <- grep("^--methods=", args, value = TRUE)
  n_rep_arg <- grep("^--n-rep=", args, value = TRUE)
  rep_start_arg <- grep("^--rep-start=", args, value = TRUE)
  full_factorial <- "--full-factorial" %in% args
  n_cores <- if (length(n_cores_arg) == 1L) {
    as.integer(sub("^--cores=", "", n_cores_arg))
  } else {
    1L
  }
  bootstrap_B <- if (length(bootstrap_arg) == 1L) {
    as.integer(sub("^--bootstrap-B=", "", bootstrap_arg))
  } else {
    NULL
  }
  method_names <- if (length(methods_arg) == 1L) {
    strsplit(sub("^--methods=", "", methods_arg), ",", fixed = TRUE)[[1L]]
  } else {
    c("UNADJ", "LADJ")
  }
  n_rep <- if (length(n_rep_arg) == 1L) {
    as.integer(sub("^--n-rep=", "", n_rep_arg))
  } else {
    NULL
  }
  rep_start <- if (length(rep_start_arg) == 1L) {
    as.integer(sub("^--rep-start=", "", rep_start_arg))
  } else {
    0L
  }

  sim <- run_simulation(
    run_full = run_full,
    n_cores = n_cores,
    scenario_design = if (full_factorial) "full_factorial" else "one_factor",
    bootstrap_B = bootstrap_B,
    method_names = method_names,
    n_rep = n_rep,
    rep_start = rep_start,
    write_outputs = write_outputs
  )
  print(sim$summary)
}
