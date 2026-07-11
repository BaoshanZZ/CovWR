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
source(file.path(.wrcov_r_dir, "config.R"))
check_sim_packages()
source(file.path(.wrcov_r_dir, "summarize_results.R"))

args <- commandArgs(trailingOnly = TRUE)
output_dir <- parse_arg(
  args,
  "output-dir",
  file.path(dirname(.wrcov_r_dir), "simulation_output_fadj_n500")
)
parts_dir <- file.path(output_dir, "parts")

read_parts <- function(pattern) {
  files <- list.files(parts_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) {
    stop("No part files matching pattern: ", pattern, call. = FALSE)
  }
  data.table::rbindlist(lapply(files, data.table::fread), use.names = TRUE, fill = TRUE)
}

raw <- read_parts("^simulation_raw_s[0-9]+\\.csv$")
truth <- read_parts("^simulation_truth_s[0-9]+\\.csv$")
scenario_grid <- read_parts("^scenario_grid_s[0-9]+\\.csv$")

summary <- summarize_simulation_results(raw)
check_simulation_outputs(raw, summary)

data.table::setorder(raw, scenario_id, rep, method, estimand)
data.table::setorder(summary, scenario_id, estimand, method)
data.table::setorder(truth, scenario_id)
data.table::setorder(scenario_grid, scenario_id)

data.table::fwrite(raw, file.path(output_dir, "simulation_raw.csv"))
data.table::fwrite(summary, file.path(output_dir, "simulation_summary.csv"))
data.table::fwrite(truth, file.path(output_dir, "simulation_truth.csv"))
data.table::fwrite(scenario_grid, file.path(output_dir, "scenario_grid.csv"))

message("Combined results written to ", output_dir)
