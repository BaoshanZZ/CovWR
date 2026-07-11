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

parse_arg <- function(args, name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 1L) {
    sub(paste0("^--", name, "="), "", hit)
  } else {
    default
  }
}

.wrcov_r_dir <- local_r_dir()
source(file.path(.wrcov_r_dir, "run_simulation.R"))

args <- commandArgs(trailingOnly = TRUE)
sid <- as.integer(parse_arg(args, "scenario-id", Sys.getenv("SLURM_ARRAY_TASK_ID")))
n_cores <- as.integer(parse_arg(args, "cores", Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
seed <- as.integer(parse_arg(args, "seed", "20260709"))
output_dir <- parse_arg(
  args,
  "output-dir",
  file.path(dirname(.wrcov_r_dir), "simulation_output_fadj_n500")
)

if (is.na(sid)) {
  stop("Missing --scenario-id or SLURM_ARRAY_TASK_ID.", call. = FALSE)
}
if (is.na(n_cores) || n_cores < 1L) {
  n_cores <- 1L
}

grid <- make_scenario_grid(run_full = TRUE, scenario_design = "one_factor")
if (!sid %in% grid$scenario_id) {
  stop("Unknown scenario_id: ", sid, call. = FALSE)
}
scenario_grid <- grid[grid[["scenario_id"]] == sid]

parts_dir <- file.path(output_dir, "parts")
dir.create(parts_dir, recursive = TRUE, showWarnings = FALSE)

message(
  "Running scenario_id=", sid,
  " with n_cores=", n_cores,
  " output_dir=", output_dir
)

sim <- run_simulation(
  run_full = TRUE,
  n_cores = n_cores,
  seed = seed,
  scenario_grid = scenario_grid,
  scenario_design = "one_factor",
  write_outputs = FALSE
)

data.table::fwrite(
  sim$raw,
  file.path(parts_dir, sprintf("simulation_raw_s%02d.csv", sid))
)
data.table::fwrite(
  sim$summary,
  file.path(parts_dir, sprintf("simulation_summary_s%02d.csv", sid))
)
data.table::fwrite(
  sim$truth,
  file.path(parts_dir, sprintf("simulation_truth_s%02d.csv", sid))
)
data.table::fwrite(
  sim$scenario_grid,
  file.path(parts_dir, sprintf("scenario_grid_s%02d.csv", sid))
)

message("Finished scenario_id=", sid)
