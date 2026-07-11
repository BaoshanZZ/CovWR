# CovWR

Simulation code for covariate-adjusted win measures with hierarchical endpoints.

The current simulation compares:

- `UNADJ`: unadjusted two-sample win-measure estimator
- `LADJ`: direct linear augmentation using raw covariates
- `FADJ`: flexible direct augmentation using a discrete Super Learner library

The flexible library contains four working learners:

- zero learner
- linear main effects
- linear main effects plus interactions
- spline terms for continuous covariates

## Repository Structure

```text
R/
  config.R                  Scenario grid and DGM settings
  dgm.R                     Data-generating mechanism
  win_kernels.R             Hierarchical win/loss kernels
  estimators.R              UNADJ, LADJ, and FADJ estimators
  truth.R                   Monte Carlo truth calculation
  summarize_results.R       Simulation summaries
  run_simulation.R          Single-process or local parallel runner
  run_scenario_batch.R      Slurm array scenario runner
  combine_batch_results.R   Slurm result combiner

submit_fadj_array.sbatch    Slurm array submission script
combine_fadj_results.sbatch Slurm combine script
```

Large simulation outputs, Slurm logs, local manuscript drafts, and local plotting helpers are
intentionally excluded from GitHub.

## R Package Requirements

The simulation uses locally installed R packages:

- `data.table`
- `future`
- `future.apply`
- `splines`

## Smoke Run

From the repository root:

```bash
Rscript -e "source('R/run_simulation.R'); run_simulation(run_full = FALSE, n_cores = 1, write_outputs = TRUE, output_dir = 'simulation_output_smoke')"
```

This runs a small number of replications for code validation.

## Full Slurm Run

The Slurm scripts are configured for the current one-factor scenario design with baseline
`N = 1000` and sample-size comparison `N = 500`.

```bash
sbatch submit_fadj_array.sbatch
```

After the array finishes, run or submit the combine script:

```bash
sbatch combine_fadj_results.sbatch
```

The combined output directory is currently configured in the scripts as:

```text
simulation_output_fadj_n1000/
```

## Summary Metrics

`simulation_summary.csv` includes:

- `estimate_mean`: Monte Carlo mean of point estimates
- `bias`: `mean(estimate - truth)`
- `ESD`: empirical standard deviation of point estimates
- `MC_SE_mean`: Monte Carlo standard error of `estimate_mean`
- `ASE`: average analytical standard error
- `CP`: empirical coverage probability
- `reject_rate`: empirical rejection rate
- `PVR`: percent variance reduction relative to `UNADJ`
