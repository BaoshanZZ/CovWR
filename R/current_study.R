# Portable entry point for the retained five-predictor study.
load_current_study <- function(project_dir=".") {
  if(!requireNamespace("data.table",quietly=TRUE)) stop("Package data.table is required.")
  library(data.table)
  data.table::setDTthreads(1L)
  project_dir <- normalizePath(project_dir,mustWork=TRUE)
  env <- new.env(parent=globalenv())
  env$project_dir <- project_dir
  files <- c("frozen_parameters.R","current_config.R","win_core.R","current_generic_kernels.R","current_kernels.R",
    "current_dgm.R","current_estimators.R","current_population.R","current_diagnostics.R","current_simulation.R","current_bootstrap.R")
  for(file in files) sys.source(file.path(project_dir,"R",file),env)
  env$source_files <- file.path(project_dir,"R",c("current_study.R","run_current_simulation.R",files))
  env$legacy_core <- env
  env$hierarchical_pair_matrices <- env$revised_pair_matrices
  env$make_population_draws <- env$full_population_draws
  env$generate_hce_outcomes <- env$staggered_outcomes
  env
}
