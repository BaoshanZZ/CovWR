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
    return(dirname(normalizePath(files[[length(files)]])))
  }
  if (dir.exists("WRCov/R")) {
    return(normalizePath("WRCov/R"))
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
                                sbw_nonnegative) {
  set.seed(rep_seed)
  dat <- generate_trial_data(
    N = scenario$N,
    e = scenario$e,
    effect = scenario$effect,
    prognosis = scenario$prognosis,
    outcome_shape = scenario$outcome_shape,
    effect_heterogeneity = scenario$effect_heterogeneity,
    endpoint_prognosis_pattern = scenario$endpoint_prognosis_pattern,
    params = params
  )

  truth_values <- truth_vector_for_scenario(truth_row)
  est <- estimate_all_methods(
    dat = dat,
    truth_values = truth_values,
    params = params,
    sbw_nonnegative = sbw_nonnegative
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
    endpoint_prognosis_pattern = scenario$endpoint_prognosis_pattern
  )]

  data.table::setcolorder(est, c(
    "scenario_id", "rep", "N", "e", "effect", "prognosis",
    "outcome_shape", "effect_heterogeneity", "endpoint_prognosis_pattern",
    "method", "estimand", "estimate", "se", "ci_low", "ci_high", "p_value",
    "reject", "truth", "covered", "tau_w", "tau_l", "learner_w", "learner_l"
  ))
  est[]
}

run_simulation <- function(run_full = FALSE,
                           n_cores = 1L,
                           seed = 20260617,
                           scenario_grid = NULL,
                           scenario_design = c("one_factor", "full_factorial"),
                           write_outputs = FALSE,
                           output_dir = file.path(dirname(.wrcov_r_dir), "simulation_output")) {
  params <- default_dgm_params()
  settings <- simulation_settings(run_full)
  scenario_design <- match.arg(scenario_design)
  grid <- if (is.null(scenario_grid)) {
    make_scenario_grid(run_full, scenario_design = scenario_design)
  } else {
    normalize_scenario_grid(scenario_grid)
  }

  message("Estimating truth values for ", nrow(grid), " scenario(s).")
  truth <- estimate_truth_grid(
    grid = grid,
    n_pairs = settings$truth_pairs,
    chunk_size = settings$truth_chunk_size,
    seed = seed,
    params = params
  )

  tasks <- data.table::CJ(
    scenario_id = grid$scenario_id,
    rep = seq_len(settings$n_rep),
    sorted = FALSE
  )
  tasks[, rep_seed := seed + scenario_id * 100000L + rep]

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
      sbw_nonnegative = settings$sbw_nonnegative
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
  full_factorial <- "--full-factorial" %in% args
  n_cores <- if (length(n_cores_arg) == 1L) {
    as.integer(sub("^--cores=", "", n_cores_arg))
  } else {
    1L
  }

  sim <- run_simulation(
    run_full = run_full,
    n_cores = n_cores,
    scenario_design = if (full_factorial) "full_factorial" else "one_factor",
    write_outputs = write_outputs
  )
  print(sim$summary)
}
