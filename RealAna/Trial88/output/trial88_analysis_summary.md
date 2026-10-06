# Trial 88 fixed-horizon hierarchical win analysis

Generated: 2026-08-22 11:07:33 EDT

## Analysis specification

- Participants: N = 200; Nutrition education n1 = 101; Control n0 = 99.
- Known simple-randomization probability: e = 0.5; observed n1/N was not substituted.
- Treatment coding: Nutrition education = 1; Control = 0.
- Fixed horizon: 1 year.
- Hierarchy: 1-year survival > AMI-free status > repeat-revascularization-free status.
- Two deaths are tied overall because death times are unavailable.
- Positive NB and WR/WO above 1 favor Nutrition education.
- Primary inference uses the subject-level Hajek IF covariance estimator; no bootstrap is used.

## Data validation and event sparsity

All four outcomes are complete binary variables, and the dictionary defines 1 as occurrence of the adverse event.

| Outcome | Nutrition_education_events | Control_events | Total_events |
| --- | --- | --- | --- |
| YP_death_1y | 5 | 7 | 12 |
| YS_ami_1y | 5 | 6 | 11 |
| YS_repeat_pci_1y | 4 | 6 | 10 |
| YS_cabg_1y | 4 | 4 | 8 |

All nine primary covariates are baseline, pretreatment variables in the dictionary and have zero missingness.

## Pairwise hierarchy summary

| N | n1 | n0 | N_pairs | Treatment_wins | Treatment_losses | Final_ties | tau_w | tau_l | tau_tie | Decided_by_death | Decided_by_death_percent | Decided_by_AMI | Decided_by_AMI_percent | Decided_by_revascularization | Decided_by_revascularization_percent | Final_ties_percent | Final_ties_at_death | Final_ties_at_revascularization |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 200.0000 | 101.0000 | 99.0000 | 9999.0000 | 1230.0000 | 983.0000 | 7786.0000 | 0.1230 | 0.0983 | 0.7787 | 1132.0000 | 11.3211 | 636.0000 | 6.3606 | 445.0000 | 4.4504 | 77.8678 | 35.0000 | 7751.0000 |

Final tie probability = 0.7787. Both elementary probabilities exceed 0.05; LWR is defined, although sparse events still warrant caution.

## Primary UNADJ and LADJ results

| Method | Measure | Estimate_log_scale | Estimate_report_scale | SE | CI_lower | CI_upper | P_value | Estimated_variance | Variance_reduction_percent | Relative_efficiency | tau_w | tau_l | tau_tie |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| UNADJ | NB |  | 0.02470 | 0.04603 | -0.06551 | 0.11491 | 0.59147 | 0.00212 | 0.00000 | 1.00000 | 0.12301 | 0.09831 | 0.77868 |
| UNADJ | WR | 0.22416 | 1.25127 | 0.41875 | 0.55070 | 2.84309 | 0.59243 | 0.17535 | 0.00000 | 1.00000 | 0.12301 | 0.09831 | 0.77868 |
| UNADJ | WO | 0.04941 | 1.05066 | 0.09211 | 0.87712 | 1.25853 | 0.59162 | 0.00848 | 0.00000 | 1.00000 | 0.12301 | 0.09831 | 0.77868 |
| LADJ | NB |  | 0.03224 | 0.04554 | -0.05701 | 0.12149 | 0.47896 | 0.00207 | 2.11739 | 1.02163 | 0.13358 | 0.10134 | 0.76509 |
| LADJ | WR | 0.27621 | 1.31813 | 0.39350 | 0.60955 | 2.85039 | 0.48272 | 0.15484 | 11.69591 | 1.13245 | 0.13358 | 0.10134 | 0.76509 |
| LADJ | WO | 0.06450 | 1.06663 | 0.09117 | 0.89209 | 1.27530 | 0.47926 | 0.00831 | 2.03328 | 1.02075 | 0.13358 | 0.10134 | 0.76509 |

For WR and WO, SE and p-values are on the log scale; reported estimates and confidence limits are exponentiated. Final adjusted measures are exact transformations of adjusted tau_w and tau_l.

## Primary efficiency comparison

| Method | Measure | Unadjusted_variance | Adjusted_variance | Variance_reduction_percent | Relative_efficiency |
| --- | --- | --- | --- | --- | --- |
| LADJ | NB | 0.00212 | 0.00207 | 2.11739 | 1.02163 |
| LADJ | WR | 0.17535 | 0.15484 | 11.69591 | 1.13245 |
| LADJ | WO | 0.00848 | 0.00831 | 2.03328 | 1.02075 |

Variance reduction is reported as computed and is not truncated at zero.

## Primary linear projection strength

| Score | Interpretation | Fixed_gradient_projection_reduction_percent |
| --- | --- | --- |
| psi_w | Elementary win IF | 4.4378 |
| psi_l | Elementary loss IF | 1.5198 |
| psi_w_plus_psi_l | Non-tie/activity direction | 4.5709 |
| psi_w_minus_psi_l | NB direction | 2.1174 |
| NB_IF_at_unadjusted_gradient | NB score using one common fixed gradient | 2.1174 |
| WR_IF_at_unadjusted_gradient | WR score using one common fixed gradient | 1.8404 |
| WO_IF_at_unadjusted_gradient | WO score using one common fixed gradient | 2.1174 |

These fixed-gradient reductions isolate how much of each IF score is removed by the linear projection. They differ in finite samples from comparing each estimator's final plug-in variance because WR gradients are evaluated at different unadjusted and adjusted elementary estimates.

## Primary baseline balance

| Covariate | Treatment_mean | Control_mean | Raw_mean_difference | Pooled_within_arm_SD | Overall_SD | Standardized_mean_difference |
| --- | --- | --- | --- | --- | --- | --- |
| X_age_years | 60.3861 | 59.5556 | 0.8306 | 9.7934 | 9.7776 | 0.0848 |
| X_sex_code | 0.6931 | 0.7071 | -0.0140 | 0.4605 | 0.4594 | -0.0304 |
| X_hypertension | 0.8218 | 0.8889 | -0.0671 | 0.3523 | 0.3530 | -0.1905 |
| X_dyslipidemia | 0.6931 | 0.8283 | -0.1352 | 0.4238 | 0.4282 | -0.3190 |
| X_diabetes | 0.4059 | 0.3636 | 0.0423 | 0.4886 | 0.4878 | 0.0866 |
| X_heart_failure | 0.4455 | 0.3434 | 0.1021 | 0.4886 | 0.4901 | 0.2090 |
| X_prior_ami | 0.4356 | 0.4242 | 0.0114 | 0.4975 | 0.4963 | 0.0229 |
| X_prior_pci | 0.1881 | 0.2424 | -0.0543 | 0.4120 | 0.4119 | -0.1318 |
| X_prior_cabg | 0.0891 | 0.1313 | -0.0422 | 0.3138 | 0.3137 | -0.1345 |

The standardized mean difference uses the pooled within-arm SD.

## Primary LADJ coefficients

| Method | Covariate | beta_w | beta_l |
| --- | --- | --- | --- |
| LADJ | X_age_years | -0.0021984 | 0.0034296 |
| LADJ | X_sex_code | 0.0909609 | -0.0452663 |
| LADJ | X_hypertension | 0.2797096 | 0.0489365 |
| LADJ | X_dyslipidemia | -0.1185233 | 0.0827201 |
| LADJ | X_diabetes | 0.0261768 | 0.0978270 |
| LADJ | X_heart_failure | -0.1668962 | 0.0248237 |
| LADJ | X_prior_ami | -0.1508699 | 0.0027418 |
| LADJ | X_prior_pci | 0.0691478 | 0.0544296 |
| LADJ | X_prior_cabg | 0.3556152 | 0.1155523 |

## Primary LADJ diagnostics

| Method | Number_features | Gram_rank | Smallest_eigenvalue | Largest_eigenvalue | Condition_number | Scaled_condition_number | Max_centering_identity_error | Max_tau_formula_error | Max_scaled_fitted_value_error | Min_adjusted_covariance_eigenvalue |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LADJ | 9.0000000 | 9.0000000 | 3.8655712 | 4757.5745717 | 1230.7558871 | 3.0357902 | 0.0000000 | 0.0000000 | 0.0000000 | 0.1392140 |

## All-covariate sensitivity coding

| Candidate_coding | Number_columns | Gram_rank | Smallest_eigenvalue | Largest_eigenvalue | Condition_number | Exact_relation | Resolution | Max_raw_column_residual_from_final_coded_span |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| All dictionary X columns; smoking dummy-coded; procedure indicators unreparameterized | 18.00000 | 17.00000 | -0.00000 | 4824.25719 |    Inf | X_balloon_angioplasty + X_stent + X_pharmacologic_stent = 1; centered columns sum to zero | Reference-code the three-level procedure category with balloon angioplasty as reference; this preserves its full two-dimensional centered span | 0.00000 |

Smoking was dummy-coded with 0 as reference. The three mutually exclusive procedure indicators were treated as one categorical variable; X_balloon_angioplasty was the reference. This removes only a redundant centered column and preserves all procedure-category information.

## All-covariate sensitivity results

| Method | Measure | Estimate_log_scale | Estimate_report_scale | SE | CI_lower | CI_upper | P_value | Variance_reduction_percent | Relative_efficiency | tau_w | tau_l | tau_tie |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LADJ sensitivity using all eligible complete baseline covariates | NB |  | 0.03975 | 0.04494 | -0.04833 | 0.12783 | 0.37644 | 4.65738 | 1.04885 | 0.13822 | 0.09847 | 0.76331 |
| LADJ sensitivity using all eligible complete baseline covariates | WR | 0.33908 | 1.40366 | 0.38999 | 0.65359 | 3.01454 | 0.38459 | 13.26408 | 1.15292 | 0.13822 | 0.09847 | 0.76331 |
| LADJ sensitivity using all eligible complete baseline covariates | WO | 0.07954 | 1.08279 | 0.09002 | 0.90764 | 1.29173 | 0.37694 | 4.47208 | 1.04681 | 0.13822 | 0.09847 | 0.76331 |

| Method | Number_features | Gram_rank | Smallest_eigenvalue | Largest_eigenvalue | Condition_number | Scaled_condition_number | Max_centering_identity_error | Max_tau_formula_error | Max_scaled_fitted_value_error | Min_adjusted_covariance_eigenvalue |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LADJ sensitivity using all eligible complete baseline covariates | 17.0000000 | 17.0000000 | 0.7928393 | 4824.2031343 | 6084.7175034 | 18.5123828 | 0.0000000 | 0.0000000 | 0.0000000 | 0.1340147 |

The final SD-scaled sensitivity Gram matrix is acceptably conditioned.

## Quality control

All 24 prespecified quality-control checks passed.

| Check | Observed | Tolerance | Passed | Details |
| --- | --- | --- | --- | --- |
| Number of pairs equals n1*n0 | 9999.00000000 | 9999.00000000 | TRUE |  |
| Each pair has win+loss+tie=1 | 0.00000000 | 0.00000000 | TRUE |  |
| Pair probabilities sum to one | 0.00000000 | 0.00000000 | TRUE |  |
| Placement means are zero | 0.00000000 | 0.00000000 | TRUE |  |
| Mean elementary IFs are zero | 0.00000000 | 0.00000000 | TRUE |  |
| Known e is exactly 0.5 | 0.50000000 | 0.50000000 | TRUE |  |
| No intercept in LADJ | 0.00000000 | 0.00000000 | TRUE |  |
| Primary LADJ uses exact prespecified covariates | 9.00000000 | 9.00000000 | TRUE |  |
| Primary adjusted tau formulas agree | 0.00000000 | 0.00000000 | TRUE |  |
| Sensitivity adjusted tau formulas agree | 0.00000000 | 0.00000000 | TRUE |  |
| Primary Gram matrix is full rank | 9.00000000 | 9.00000000 | TRUE |  |
| Sensitivity final coding is full rank | 17.00000000 | 17.00000000 | TRUE |  |
| Raw all-X rank deficiency was detected | 17.00000000 | 17.00000000 | TRUE | X_balloon_angioplasty + X_stent + X_pharmacologic_stent = 1; centered columns sum to zero |
| Final sensitivity coding preserves raw all-X centered span | 0.00000000 | 0.00000000 | TRUE |  |
| Primary covariance is symmetric | 0.00000000 | 0.00000000 | TRUE |  |
| Sensitivity covariance is symmetric | 0.00000000 | 0.00000000 | TRUE |  |
| Primary covariance is positive semidefinite | 0.13921396 | -0.00000000 | TRUE |  |
| Sensitivity covariance is positive semidefinite | 0.13401472 | -0.00000000 | TRUE |  |
| All four outcomes are complete | 0.00000000 | 0.00000000 | TRUE |  |
| All four outcomes use 0/1 adverse-event coding | 1.00000000 | 1.00000000 | TRUE |  |
| Adverse components reproduce composite event | 0.00000000 | 0.00000000 | TRUE |  |
| All randomized subjects retained | 200.00000000 | 200.00000000 | TRUE |  |
| tau_w is positive for LWR | 0.12301230 | 0.00000000 | TRUE |  |
| tau_l is positive for LWR | 0.09830983 | 0.00000000 | TRUE |  |

## Interpretation cautions

- Cardiovascular events are sparse and most treatment-control pairs are expected to remain tied; the effective information is determined by subjects, not by treating all pairs as independent.
- Wald inference relies on the first-order subject-level influence-function approximation.
- The all-covariate analysis is a sensitivity analysis, not the primary analysis; it uses many predictors relative to N=200.
- The observed variance reduction may be negative even though the projection no-loss statement is asymptotic.
