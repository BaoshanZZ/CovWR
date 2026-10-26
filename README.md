# CovWR

Reproducible R code for unadjusted and linearly covariate-adjusted marginal win
statistics with fixed-horizon hierarchical composite endpoints: net benefit (NB),
log win ratio (LWR), and log win odds (LWO).

This release uses the retained **five-predictor simulation**. It replaces the older
N=1000/N=500, 20,000-replication one-factor-at-a-time study previously stored here.
The N=200 bootstrap comparison is a separate variance-estimation sensitivity
analysis, not a replacement for the primary DGP.

## Contents

```text
R/
  current_study.R              Portable loader
  run_current_simulation.R     Smoke/full simulation, bootstrap and population diagnostics
  validate_current_study.R     Kernel, IF, projection, transformation and bootstrap checks
  frozen_parameters.R         Exact retained model parameters; no local RDS dependency
  current_*.R, win_core.R      Reused statistical implementation
simulation_results/
  primary/                    Final primary aggregate results and population benchmarks
  bootstrap/                  N=200 paired bootstrap-versus-IF summaries
  code_provenance.csv         Functions extracted from the retained implementation
figures/                      Six final EPS figures only
RealAna/Trial3/                AKI alert trial code, dictionary and aggregate outputs
```

Raw clinical datasets, MC replication outputs, checkpoints, scheduler logs,
batch/bash files, manuscript TeX, and plotting scripts are not published.

## Statistical Method

The methods first estimate elementary treatment win/loss probabilities using
two-sample U-statistics and centered subject-level placements. LADJ projects the
two elementary IFs onto `(A-e) * (X_W-Xbar_W)` with known `e=0.5`, no intercept,
and unregularized QR least squares. The adjusted elementary vector is transformed
exactly to NB, LWR and LWO. SEs use the empirical elementary-IF covariance and the
delta method at the corresponding estimated win/loss vector.

The code checks rank and numerical conditioning; it does not silently regularize
or clip probabilities. No propensity estimation, SBW, entropy balancing, or
flexible learner is part of the retained simulation.

## Primary Design

- Sample sizes: **N=200,600,800,1200**; simple Bernoulli randomization, `e=0.5`.
- **10,000 null and 10,000 alternative trials per N**: 80,000 trials total.
- Common horizon: **S=28 days**, with no dropout before S.
- Five correlated prognostic variables from a latent Gaussian vector with
  off-diagonal correlation 0.2. X1/X3/X5 are continuous; X2/X4 are binary with
  prevalences 0.50/0.35. Observed correlations after dichotomization differ from 0.2.
- D1: terminal Weibull time to event; D2: rounded/bounded Gaussian duration;
  D3: positive lognormal measurement, structurally unavailable after the terminal event.
- Longer survival and lower D2/D3 are favorable. Among two survivors, D2 uses
  a two-unit tie margin and D3 a larger-to-smaller ratio margin of 1.25.
- Partially overlapping direct predictors: X1-X4 for D1, X2-X4 for D2, X3-X5 for D3.
- Null treatment coefficients are zero. The fixed alternative vector is
  `1.3125 * c(log(0.85), -0.50, log(0.95))` for every N and working set.
- X7 and X9 are independent Bernoulli noise with prevalences 0.5 and 0.3;
  X6/X8/X10-X15 are independent standard normal noise. All noise coefficients are zero.

### Working Sets

| Code | Label | Covariates | p |
|---|---|---|---:|
| W0 | Unadjusted | None | 0 |
| W1 | Shared | X3, X4 | 2 |
| W2 | D1-pair, exploratory | X1, X2 | 2 |
| W3 | D1-prog | X1-X4 | 4 |
| W4 | D2-prog | X2-X4 | 3 |
| W5 | D3-prog | X3-X5 | 3 |
| W6 | Union | X1-X5 | 5 |
| W7 | Union+5 | X1-X10 | 10 |
| W8 | Union+10 | X1-X15 | 15 |

All nine sets use the same generated trial, kernels and elementary IFs within a
replication. The manuscript/figures report seven adjusted sets plus unadjusted;
the exploratory W2 is retained in numerical CSVs but not displayed. `Union` was
previously labeled `Full`; `Dq-prog` was labeled `Dq-complete`.

## Bootstrap Sensitivity Analysis

**Pooled subject bootstrap:** uses the original primary N=200 trials, 10,000 per
null/alternative setting, with K=2000 resamples per trial. Whole subjects are
resampled with treatment labels retained; arm sizes may vary while e remains 0.5.
Placements, elementary IFs, centering and adjustment coefficients are refitted.
The SD of bootstrap point estimates replaces the analytical SE in the same
normal Wald procedure; original point estimates remain unchanged. This is not
percentile, BCa, or bootstrap-t inference. All reported bootstrap samples succeeded.

## Reproduce

Run from the repository root with R 4.4 or later. `data.table` is required;
`survival` is additionally required for component prediction diagnostics.
No package installation is performed by these scripts.

```bash
Rscript R/validate_current_study.R
Rscript R/run_current_simulation.R --mode=smoke
Rscript R/run_current_simulation.R --stage=bootstrap --mode=smoke
```

Smoke mode uses five trials per requested cell and K=20 for bootstrap. Its outputs
are execution checks, not precise operating-characteristic estimates.

```bash
Rscript R/run_current_simulation.R --experiment=primary --mode=full --cores=4
Rscript R/run_current_simulation.R --stage=bootstrap --mode=full --effects=null,alternative
Rscript R/run_current_simulation.R --stage=population --experiment=primary --mode=full
```

The full bootstrap run is expensive. Smaller runs can use `--reps=100`, `--K=200`,
`--samples=200,800`, or `--effects=null` where applicable. Multicore simulation
uses Unix fork workers; use `--cores=1` on Windows. Bootstrap runs one trial at a time.
Outputs default to ignored `simulation_runs/`; `--output=...` can select another
local destination. Checkpoints are resumable only when settings/streams match.
The saved `simulation_results/` summaries are protected from overwriting.

Compact random-stream schedules retain the first exact L'Ecuyer-CMRG state of each
10,000-trial block. Trials use the retained master seed 2026102306; bootstrap uses
2026102503. Different N are independently generated, not nested subsamples.
Changing the requested number of trials does not change the first replicates.
Retained parameter and stream files, rather than a newly selected seed, define
reproduction. Exact last-bit agreement can depend on R and numerical-library versions.

## Results and Benchmarks

The primary `simulation_summary.csv` contains 216 rows covering
N x effect x W x NB/LWR/LWO. It reports bias, ESD, average analytical SE (ASE),
ASE/ESD, empirical and estimated variance reduction, coverage, rejection rates,
successful/failed counts, and pointwise 95% Wilson MC intervals. Empirical PVR
compares ESD squared with unadjusted within the same N/effect/DGP.

`simulation_truth.csv` records independent-pair win/loss/tie probabilities, NB/LWR/LWO,
truth MCSEs and component-resolution probabilities. Alternative truths use 10^7
independent pairs; null targets are exactly zero by symmetry. Population prediction
uses separate 100,000-subject training/validation controls. Population IF-projection
PVR uses 100,000 training/evaluation subjects and independent 100,000-subject
reference banks per arm/split; it is a linear-class benchmark, not outcome-regression
performance or an unrestricted efficiency bound. No new parameter calibration is
performed by the population stage; rerunning finite-MC benchmarks can give slightly
different approximations from retained ones.

Finite-sample analytical SE underestimation can occur, especially at N=200 with
noise covariates. Higher unadjusted-for-size power is not necessarily greater
efficiency. No-loss is an asymptotic projection property, not a guarantee that each
finite-sample variance estimate or empirical PVR improves. These are retained
post-exploration configurations, not prospectively registered simulation designs;
MC intervals do not adjust for earlier design or seed comparisons.

## AKI Alert Trial Application

`RealAna/Trial3/` corresponds to the electronic health record alert trial reported
by Wilson et al., [BMJ 2021;372:m4786](https://doi.org/10.1136/bmj.m4786).
The analysis uses N=6030 randomized participants, alert versus usual care,
known e=0.5 and binary 14-day death > dialysis > AKI-progression status.
It compares unadjusted and LADJ inference without primary bootstrap inference.

The de-identified RDS is not redistributed. Place the source data at
`RealAna/Trial3/trial3.rds`, beside the supplied dictionary, then run:

```bash
Rscript RealAna/Trial3/analysis_aki_alert_win_ladj.R
Rscript RealAna/Trial3/build_baseline_table.R
```

All three analyzed outcomes and seven selected baseline variables are complete.
Continuous age (507 missing) and creatinine (two missing) were excluded as variables,
not used to remove participants; the complete age-over-90 indicator was used.
No imputation was performed. The dictionary does not map sex codes 0/1 to male/female;
baseline summaries retain `Sex (code 1)` pending confirmation. Public outputs are
aggregate CSV/Markdown summaries and QC, not participant-level records.
