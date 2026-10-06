# CovWR

R implementation of unadjusted and linearly covariate-adjusted win measures for
fixed-horizon hierarchical composite endpoints.

The primary implementation estimates the elementary treatment win and loss probabilities,
then transforms them to:

- net benefit (NB);
- log win ratio (LWR) and win ratio (WR);
- log win odds (LWO) and win odds (WO).

The proposed linear adjustment (LADJ) directly projects the elementary influence
functions on sample-centered baseline covariates using the known randomization probability.
The final adjusted win measures are exact transformations of the adjusted win-loss vector.

## Repository Structure

~~~text
R/
  config.R
  dgm.R
  win_kernels.R
  estimators.R
  truth.R
  summarize_results.R
  run_simulation.R
  generate_final_lwr_figures.R
  generate_simulation_results_tables.R
  validate_method_formulas.R

simulation_results_n1000_mc20000/
  scenario_grid.csv
  simulation_summary.csv
  simulation_truth.csv
  main_report/
  supplementary_report/

RealAna/
  Trial3/
    analysis_aki_alert_win_ladj.R
    trial3_dictionary.csv
    output/
~~~

Cluster submission scripts, intermediate simulation chunks, and scheduler logs are not
versioned.

## Primary Simulation

The final simulation uses a one-factor-at-a-time design with:

- reference sample size N = 1000;
- sample-size comparison N = 500;
- balanced randomization except for the 1:2 allocation scenario;
- 20,000 Monte Carlo replications per scenario;
- 10,000,000 independently generated treatment-control pairs for alternative truth;
- UNADJ and LADJ as the primary methods.

The final aggregate results and manuscript-ready figures/tables are stored in
simulation_results_n1000_mc20000/. Replication-level raw results and scenario chunk files
are intentionally omitted.

### Smoke test

~~~bash
Rscript -e "source('R/run_simulation.R'); run_simulation(
  run_full = FALSE,
  n_cores = 1,
  write_outputs = TRUE,
  output_dir = 'simulation_output_smoke'
)"
~~~

### Regenerate final figures

~~~bash
Rscript R/generate_final_lwr_figures.R \
  --input-dir=simulation_results_n1000_mc20000 \
  --output-dir=simulation_results_n1000_mc20000/main_report
~~~

### Regenerate result tables

~~~bash
Rscript R/generate_simulation_results_tables.R \
  --input-dir=simulation_results_n1000_mc20000 \
  --output-dir=simulation_results_n1000_mc20000/main_report \
  --estimands=LWR

Rscript R/generate_simulation_results_tables.R \
  --input-dir=simulation_results_n1000_mc20000 \
  --output-dir=simulation_results_n1000_mc20000/supplementary_report \
  --estimands=NB,LWO
~~~

## Real-Data Analyses

The real-data directory contains analysis scripts, data dictionaries, numerical outputs,
quality-control results, and manuscript-ready material.

The AKI alert trial analysis is run with:

~~~bash
Rscript RealAna/Trial3/analysis_aki_alert_win_ladj.R
~~~

The .rds source datasets are not distributed through this repository. To reproduce an
analysis, place the corresponding RDS file beside its dictionary using the filename
expected by the script.

Trial 3 corresponds to the multicenter randomized clinical trial of electronic health
record alerts for acute kidney injury reported by Wilson et al. (BMJ 2021;372:m4786;
doi:10.1136/bmj.m4786).

## Validation

Formula-level checks can be run with:

~~~bash
Rscript R/validate_method_formulas.R
~~~

The script verifies the elementary influence-function construction, closed-form LADJ
solution, exact win-measure transformations, and covariance delta-method calculations.

## Package Requirements

Core simulation:

- data.table
- future
- future.apply
- splines

Figures:

- ggplot2
- scales
- ragg

The real-data scripts use base R and do not require bootstrap inference.
