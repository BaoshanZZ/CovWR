# Trial 23 fixed-horizon hierarchical win analysis

Generated: 2026-08-22 10:31:56 EDT

## Analysis specification

- Participants: N = 92; rhTM n1 = 47; No rhTM n0 = 45.
- Known simple-randomization probability: e = 0.5.
- Treatment coding: rhTM = 1; No rhTM = 0.
- YP_event_28d value denoting death/event: `1`.
- Hierarchy: 28-day survival > Day-3 change in SOFA > Day-3 change in D-dimer.
- Lower Day-3 delta SOFA and lower Day-3 delta D-dimer are favorable.
- Positive NB and WR/WO above 1 favor rhTM.
- Primary inference uses the subject-level Hajek IF covariance estimator; no bootstrap is used.

## Data validation

- Complete Day-28 event/time records: yes (92/92).
- Day-28 survivors: 75; all have Day-3 SOFA and D-dimer.
- Day-28 deaths: 17. Missing Day-3 SOFA = 6; missing Day-3 D-dimer = 7.
- Day-28 deaths by arm: rhTM = 8; No rhTM = 9.
- No hierarchical outcome was imputed and no randomized subject was removed.
- Primary LADJ covariates: X_SOFA_0d, X_SIRS_0d, X_Ddimer_0d, X_platelet_0d.
- Optional X_TM sensitivity: pooled median = 5.750; missing indicator count = 6.

## Pairwise hierarchy summary

| HCE_definition | N | n1 | n0 | N_pairs | Treatment_wins | Treatment_losses | Final_ties | tau_w | tau_l | tau_tie | Decided_by_survival | Decided_by_survival_percent | Decided_by_SOFA | Decided_by_SOFA_percent | Decided_by_Ddimer | Decided_by_Ddimer_percent | Final_ties_percent | Final_ties_at_survival | Final_ties_at_Ddimer | Reached_Ddimer_tier | SOFA_pairs_reassigned_to_Ddimer_vs_threshold_free | SOFA_pairs_reassigned_percent_of_all_pairs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Threshold-free SOFA | 92.0000 | 47.0000 | 45.0000 | 2115.0000 | 1147.0000 | 963.0000 | 5.0000 | 0.5423 | 0.4553 | 0.0024 | 706.0000 | 33.3806 | 1207.0000 | 57.0686 | 197.0000 | 9.3144 | 0.2364 | 5.0000 | 0.0000 | 197.0000 | 0.0000 | 0.0000 |
| 2-point SOFA margin | 92.0000 | 47.0000 | 45.0000 | 2115.0000 | 1155.0000 | 955.0000 | 5.0000 | 0.5461 | 0.4515 | 0.0024 | 706.0000 | 33.3806 | 838.0000 | 39.6217 | 566.0000 | 26.7612 | 0.2364 | 5.0000 | 0.0000 | 566.0000 | 369.0000 | 17.4468 |

Under the 2-point SOFA margin, 369 pairs that would have been decided by threshold-free SOFA were passed to D-dimer (17.45% of all pairs).

## Primary UNADJ and LADJ results

| HCE_definition | Method | Measure | Estimate_log_scale | Estimate_report_scale | SE | CI_lower | CI_upper | P_value | Estimated_variance | Variance_reduction_percent | Relative_efficiency |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Threshold-free SOFA | UNADJ | NB |  | 0.0870 | 0.1212 | -0.1505 | 0.3245 | 0.4728 | 0.0147 | 0.0000 | 1.0000 |
| Threshold-free SOFA | UNADJ | WR | 0.1749 | 1.1911 | 0.2448 | 0.7372 | 1.9244 | 0.4750 | 0.0599 | 0.0000 | 1.0000 |
| Threshold-free SOFA | UNADJ | WO | 0.1744 | 1.1906 | 0.2442 | 0.7377 | 1.9214 | 0.4750 | 0.0596 | 0.0000 | 1.0000 |
| Threshold-free SOFA | LADJ | NB |  | 0.0786 | 0.1206 | -0.1577 | 0.3149 | 0.5145 | 0.0145 | 0.9763 | 1.0099 |
| Threshold-free SOFA | LADJ | WR | 0.1579 | 1.1710 | 0.2432 | 0.7271 | 1.8860 | 0.5162 | 0.0591 | 1.3061 | 1.0132 |
| Threshold-free SOFA | LADJ | WO | 0.1575 | 1.1706 | 0.2427 | 0.7275 | 1.8835 | 0.5162 | 0.0589 | 1.2532 | 1.0127 |
| 2-point SOFA margin | UNADJ | NB |  | 0.0946 | 0.1198 | -0.1401 | 0.3293 | 0.4297 | 0.0143 | 0.0000 | 1.0000 |
| 2-point SOFA margin | UNADJ | WR | 0.1901 | 1.2094 | 0.2422 | 0.7523 | 1.9443 | 0.4325 | 0.0587 | 0.0000 | 1.0000 |
| 2-point SOFA margin | UNADJ | WO | 0.1897 | 1.2089 | 0.2417 | 0.7528 | 1.9413 | 0.4325 | 0.0584 | 0.0000 | 1.0000 |
| 2-point SOFA margin | LADJ | NB |  | 0.0872 | 0.1195 | -0.1470 | 0.3213 | 0.4656 | 0.0143 | 0.4826 | 1.0048 |
| 2-point SOFA margin | LADJ | WR | 0.1751 | 1.1914 | 0.2413 | 0.7425 | 1.9117 | 0.4679 | 0.0582 | 0.8034 | 1.0081 |
| 2-point SOFA margin | LADJ | WO | 0.1748 | 1.1910 | 0.2408 | 0.7430 | 1.9091 | 0.4679 | 0.0580 | 0.7520 | 1.0076 |

For WR and WO, SE and p-values are on the log scale; estimates and displayed confidence limits are exponentiated. For NB, inference is on the NB scale.

## Primary efficiency comparison

| HCE_definition | Method | Measure | Unadjusted_variance | Adjusted_variance | Variance_reduction_percent | Relative_efficiency |
| --- | --- | --- | --- | --- | --- | --- |
| Threshold-free SOFA | LADJ | NB | 0.0147 | 0.0145 | 0.9763 | 1.0099 |
| Threshold-free SOFA | LADJ | WR | 0.0599 | 0.0591 | 1.3061 | 1.0132 |
| Threshold-free SOFA | LADJ | WO | 0.0596 | 0.0589 | 1.2532 | 1.0127 |
| 2-point SOFA margin | LADJ | NB | 0.0143 | 0.0143 | 0.4826 | 1.0048 |
| 2-point SOFA margin | LADJ | WR | 0.0587 | 0.0582 | 0.8034 | 1.0081 |
| 2-point SOFA margin | LADJ | WO | 0.0584 | 0.0580 | 0.7520 | 1.0076 |

Variance reduction is reported as computed and is not truncated at zero.

## Baseline balance

| Covariate | Treatment_mean | Control_mean | Raw_mean_difference | Pooled_within_arm_SD | Overall_SD | Standardized_mean_difference |
| --- | --- | --- | --- | --- | --- | --- |
| X_SOFA_0d | 7.3404 | 8.0889 | -0.7485 | 3.9584 | 3.9545 | -0.1891 |
| X_SIRS_0d | 3.0638 | 3.2222 | -0.1584 | 0.8068 | 0.8063 | -0.1963 |
| X_Ddimer_0d | 30.1532 | 39.4689 | -9.3157 | 44.6891 | 44.6888 | -0.2085 |
| X_platelet_0d | 13.4489 | 13.3422 | 0.1067 | 7.4931 | 7.4520 | 0.0142 |

The standardized mean difference uses the pooled within-arm SD. Coefficients describe precision augmentation and are not causal effects.

## Primary LADJ coefficients

| HCE_definition | Method | Covariate | beta_w | beta_l |
| --- | --- | --- | --- | --- |
| Threshold-free SOFA | LADJ | X_SOFA_0d | 0.015002 | -0.015465 |
| Threshold-free SOFA | LADJ | X_SIRS_0d | -0.110256 | 0.112175 |
| Threshold-free SOFA | LADJ | X_Ddimer_0d | -0.001011 | 0.001122 |
| Threshold-free SOFA | LADJ | X_platelet_0d | 0.006065 | -0.006111 |
| 2-point SOFA margin | LADJ | X_SOFA_0d | -0.002081 | 0.001618 |
| 2-point SOFA margin | LADJ | X_SIRS_0d | -0.069590 | 0.071509 |
| 2-point SOFA margin | LADJ | X_Ddimer_0d | -0.000088 | 0.000199 |
| 2-point SOFA margin | LADJ | X_platelet_0d | 0.008498 | -0.008545 |

## LADJ numerical diagnostics

| HCE_definition | Method | Number_features | Gram_rank | Smallest_eigenvalue | Largest_eigenvalue | Condition_number | Scaled_condition_number | Max_centering_identity_error | Max_tau_formula_error | Max_scaled_fitted_value_error |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Threshold-free SOFA | LADJ | 4.000000 | 4.000000 | 14.274643 | 45442.224857 | 3183.422861 | 1.862065 | 0.000000 | 0.000000 | 0.000000 |
| 2-point SOFA margin | LADJ | 4.000000 | 4.000000 | 14.274643 | 45442.224857 | 3183.422861 | 1.862065 | 0.000000 | 0.000000 | 0.000000 |

Covariate rescaling to pooled SDs produced the same fitted augmentation values up to the reported numerical tolerance.

## Optional baseline-TM sensitivity

| HCE_definition | Method | Measure | Estimate_log_scale | Estimate_report_scale | SE | CI_lower | CI_upper | P_value | Variance_reduction_percent | Relative_efficiency |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Threshold-free SOFA | LADJ sensitivity including baseline TM with pooled imputation + missing indicator | NB |  | 0.1259 | 0.1172 | -0.1038 | 0.3556 | 0.2826 | 6.4762 | 1.0692 |
| Threshold-free SOFA | LADJ sensitivity including baseline TM with pooled imputation + missing indicator | WR | 0.2538 | 1.2889 | 0.2387 | 0.8073 | 2.0578 | 0.2877 | 4.8957 | 1.0515 |
| Threshold-free SOFA | LADJ sensitivity including baseline TM with pooled imputation + missing indicator | WO | 0.2532 | 1.2881 | 0.2382 | 0.8077 | 2.0543 | 0.2878 | 4.8948 | 1.0515 |
| 2-point SOFA margin | LADJ sensitivity including baseline TM with pooled imputation + missing indicator | NB |  | 0.1256 | 0.1168 | -0.1034 | 0.3546 | 0.2825 | 4.7908 | 1.0503 |
| 2-point SOFA margin | LADJ sensitivity including baseline TM with pooled imputation + missing indicator | WR | 0.2531 | 1.2880 | 0.2380 | 0.8079 | 2.0536 | 0.2876 | 3.4599 | 1.0358 |
| 2-point SOFA margin | LADJ sensitivity including baseline TM with pooled imputation + missing indicator | WO | 0.2525 | 1.2872 | 0.2374 | 0.8083 | 2.0501 | 0.2876 | 3.4652 | 1.0359 |

This sensitivity analysis is separate from the primary result. X_TM was imputed using the pooled observed median without treatment or outcome information, and a missingness indicator was included.

## Quality control

All 31 prespecified quality-control checks passed.

| HCE_definition | Check | Observed | Tolerance | Passed |
| --- | --- | --- | --- | --- |
| Threshold-free SOFA | Number of pairs equals n1*n0 | 2115.00000000 | 2115.00000000 | TRUE |
| Threshold-free SOFA | Each pair has win+loss+tie=1 | 0.00000000 | 0.00000000 | TRUE |
| Threshold-free SOFA | Pair probabilities sum to one | 0.00000000 | 0.00000000 | TRUE |
| Threshold-free SOFA | Placement means are zero | 0.00000000 | 0.00000000 | TRUE |
| Threshold-free SOFA | Mean elementary IFs are zero | 0.00000000 | 0.00000000 | TRUE |
| Threshold-free SOFA | Known e is exactly 0.5 | 0.50000000 | 0.50000000 | TRUE |
| Threshold-free SOFA | No intercept in primary LADJ | 0.00000000 | 0.00000000 | TRUE |
| Threshold-free SOFA | Primary LADJ uses exactly four prespecified covariates | 4.00000000 | 4.00000000 | TRUE |
| Threshold-free SOFA | Equivalent adjusted tau formulas agree | 0.00000000 | 0.00000000 | TRUE |
| Threshold-free SOFA | Primary Gram matrix is full rank | 4.00000000 | 4.00000000 | TRUE |
| Threshold-free SOFA | Adjusted covariance is symmetric | 0.00000000 | 0.00000000 | TRUE |
| Threshold-free SOFA | Adjusted covariance is positive semidefinite | 0.00014643 | -0.00000000 | TRUE |
| Threshold-free SOFA | No survivor missing Day-3 SOFA | 0.00000000 | 0.00000000 | TRUE |
| Threshold-free SOFA | No pair reaching D-dimer has missing D-dimer | 0.00000000 | 0.00000000 | TRUE |
| Threshold-free SOFA | All randomized subjects retained | 92.00000000 | 92.00000000 | TRUE |
| 2-point SOFA margin | Number of pairs equals n1*n0 | 2115.00000000 | 2115.00000000 | TRUE |
| 2-point SOFA margin | Each pair has win+loss+tie=1 | 0.00000000 | 0.00000000 | TRUE |
| 2-point SOFA margin | Pair probabilities sum to one | 0.00000000 | 0.00000000 | TRUE |
| 2-point SOFA margin | Placement means are zero | 0.00000000 | 0.00000000 | TRUE |
| 2-point SOFA margin | Mean elementary IFs are zero | 0.00000000 | 0.00000000 | TRUE |
| 2-point SOFA margin | Known e is exactly 0.5 | 0.50000000 | 0.50000000 | TRUE |
| 2-point SOFA margin | No intercept in primary LADJ | 0.00000000 | 0.00000000 | TRUE |
| 2-point SOFA margin | Primary LADJ uses exactly four prespecified covariates | 4.00000000 | 4.00000000 | TRUE |
| 2-point SOFA margin | Equivalent adjusted tau formulas agree | 0.00000000 | 0.00000000 | TRUE |
| 2-point SOFA margin | Primary Gram matrix is full rank | 4.00000000 | 4.00000000 | TRUE |
| 2-point SOFA margin | Adjusted covariance is symmetric | 0.00000000 | 0.00000000 | TRUE |
| 2-point SOFA margin | Adjusted covariance is positive semidefinite | 0.00014652 | -0.00000000 | TRUE |
| 2-point SOFA margin | No survivor missing Day-3 SOFA | 0.00000000 | 0.00000000 | TRUE |
| 2-point SOFA margin | No pair reaching D-dimer has missing D-dimer | 0.00000000 | 0.00000000 | TRUE |
| 2-point SOFA margin | All randomized subjects retained | 92.00000000 | 92.00000000 | TRUE |
| Across HCE definitions | Same subjects and baseline procedure; only SOFA rule changes | 1.00000000 | 1.00000000 | TRUE |

## Interpretation cautions

- This is a small trial; Wald inference relies on a first-order subject-level influence-function approximation.
- The 2-point SOFA rule is a pairwise clinical relevance margin, not an individual responder definition.
- Estimated variance reduction in this realized dataset need not be nonnegative even though LADJ has an asymptotic no-loss property under the stated projection conditions.
- Clinical suitability also depends on whether the prespecified hierarchy and direction of the delta outcomes are defensible for the Trial 23 scientific question.
