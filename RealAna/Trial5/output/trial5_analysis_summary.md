# Trial 5: 14-day hierarchical win analysis

## Verified mappings

| Role | Exact_variable | Coding_or_use |
| --- | --- | --- |
| Treatment | Treatment | A=1: Alert; A=0: Usual Care |
| 14-day death | YS_death_14d | 1=death occurrence; 0=no death |
| 14-day dialysis | YS_dialysis_14d | 1=dialysis occurrence; 0=no dialysis |
| 14-day AKI progression | YS_aki_progression_14d | 1=AKI progression occurrence; 0=no progression |
| Hospital / study site | X_hospital_0m | Baseline LADJ covariate |
| Age (continuous) | X_age_0m | Baseline LADJ covariate |
| Sex | X_sex_0m | Baseline LADJ covariate |
| ICU status | X_icu_0m | Baseline LADJ covariate |
| Baseline serum creatinine | X_creat_0m | Baseline LADJ covariate |
| Elixhauser comorbidity score | X_elx_score_pmhx_0m | Baseline LADJ covariate |
| Diabetes | X_diabetes_elx_pmhx_0m | Baseline LADJ covariate |
| Congestive heart failure | X_chf_elx_pmhx_0m | Baseline LADJ covariate |

## Analysis population

| Trial | Population | Role | N | Treatment_code_A | Known_randomization_probability |
| --- | --- | --- | --- | --- | --- |
| Trial 5 | All randomized | Total | 5060.0000 |  | 0.5000 |
| Trial 5 | Usual Care | Control | 2528.0000 | 0.0000 | 0.5000 |
| Trial 5 | Alert | Intervention | 2532.0000 | 1.0000 | 0.5000 |

## Outcome missingness

| Outcome | Variable | Missing |
| --- | --- | --- |
| 14-day death | YS_death_14d | 0.0000 |
| 14-day dialysis | YS_dialysis_14d | 0.0000 |
| 14-day AKI progression | YS_aki_progression_14d | 0.0000 |

## Table A: outcomes by randomized arm

| Arm | Role | N | Death_n | Death_percent | Dialysis_n | Dialysis_percent | AKI_progression_n | AKI_progression_percent | Severity_death_n | Severity_death_percent | Severity_alive_dialysis_n | Severity_alive_dialysis_percent | Severity_alive_no_dialysis_AKI_n | Severity_alive_no_dialysis_AKI_percent | Severity_none_n | Severity_none_percent |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Usual Care | Control | 2528.0000 | 282.0000 | 11.1551 | 127.0000 | 5.0237 | 505.0000 | 19.9763 | 282.0000 | 11.1551 | 71.0000 | 2.8085 | 286.0000 | 11.3133 | 1889.0000 | 74.7231 |
| Alert | Intervention | 2532.0000 | 253.0000 | 9.9921 | 123.0000 | 4.8578 | 475.0000 | 18.7599 | 253.0000 | 9.9921 | 65.0000 | 2.5671 | 267.0000 | 10.5450 | 1947.0000 | 76.8957 |

## Unadjusted pairwise comparison

| Trial | Horizon_days | N | n1 | n0 | N_pairs | Treatment_wins | Treatment_losses | Final_ties | tau_w | tau_l | tau_tie | Decided_at_death_n | Decided_at_death_percent | Decided_at_dialysis_n | Decided_at_dialysis_percent | Decided_at_AKI_progression_n | Decided_at_AKI_progression_percent | Final_ties_percent |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Trial 5 | 14.0000 | 5060.0000 | 2532.0000 | 2528.0000 | 6400896.0000 | 1377864.0000 | 1238454.0000 | 3784578.0000 | 0.2153 | 0.1935 | 0.5913 | 1210916.0000 | 18.9179 | 322597.0000 | 5.0399 | 1082805.0000 | 16.9165 | 59.1258 |

## Table B: prespecified baseline covariate audit

| Clinical_covariate | Variable | Missing | Selected | Reason |
| --- | --- | --- | --- | --- |
| Hospital / study site | X_hospital_0m | 0.0000 | TRUE | Included: verified baseline and fully observed. |
| Age (continuous) | X_age_0m | 0.0000 | TRUE | Included: continuous baseline age is complete. |
| Age category (fallback) | X_age_over_90_0m | 0.0000 | FALSE | Not included because complete continuous age was available. |
| Sex | X_sex_0m | 0.0000 | TRUE | Included: verified baseline and fully observed. |
| ICU status | X_icu_0m | 0.0000 | TRUE | Included: verified baseline and fully observed. |
| Baseline serum creatinine | X_creat_0m | 0.0000 | TRUE | Included: verified baseline and fully observed. |
| Baseline AKI stage |  |  | FALSE | No true baseline AKI-stage variable is present; YS_aki_stage_14d is post-randomization and was excluded. |
| Elixhauser comorbidity score | X_elx_score_pmhx_0m | 0.0000 | TRUE | Included: verified baseline and fully observed. |
| Diabetes | X_diabetes_elx_pmhx_0m | 0.0000 | TRUE | Included: verified baseline and fully observed. |
| Congestive heart failure | X_chf_elx_pmhx_0m | 0.0000 | TRUE | Included: verified baseline and fully observed. |

**Final LADJ raw covariates:** X_hospital_0m, X_age_0m, X_sex_0m, X_icu_0m, X_creat_0m, X_elx_score_pmhx_0m, X_diabetes_elx_pmhx_0m, X_chf_elx_pmhx_0m.

**Final LADJ design columns after reference dummy coding:** X_hospital_0m2, X_hospital_0m3, X_hospital_0m4, X_age_0m, X_sex_0m1, X_icu_0m1, X_creat_0m, X_elx_score_pmhx_0m, X_diabetes_elx_pmhx_0m1, X_chf_elx_pmhx_0m1.

## Table C: unadjusted and LADJ estimates

| Trial | Estimand | Method | Estimate | ASE | CI_lower | CI_upper | P_value | ASE_ratio_LADJ_to_UNADJ | Relative_efficiency_vs_UNADJ | Percent_variance_reduction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Trial 5 | NB | UNADJ | 0.021780 | 0.012174 | -0.002081 | 0.045641 | 0.073612 | 1.000000 | 1.000000 | 0.000000 |
| Trial 5 | NB | LADJ | 0.016692 | 0.011611 | -0.006066 | 0.039449 | 0.150558 | 0.953748 | 1.099342 | 9.036519 |
| Trial 5 | LWR | UNADJ | 0.106671 | 0.059671 | -0.010282 | 0.223624 | 0.073833 | 1.000000 | 1.000000 | 0.000000 |
| Trial 5 | LWR | LADJ | 0.081788 | 0.056910 | -0.029752 | 0.193329 | 0.150671 | 0.953722 | 1.099402 | 9.041501 |
| Trial 5 | WR | UNADJ | 1.112568 | 0.066388 | 0.989770 | 1.250600 | 0.073833 | 1.000000 | 1.000000 | 0.000000 |
| Trial 5 | WR | LADJ | 1.085226 | 0.061760 | 0.970686 | 1.213282 | 0.150671 | 0.930284 | 1.155498 | 13.457213 |
| Trial 5 | LWO | UNADJ | 0.043566 | 0.024360 | -0.004178 | 0.091311 | 0.073703 | 1.000000 | 1.000000 | 0.000000 |
| Trial 5 | LWO | LADJ | 0.033386 | 0.023229 | -0.012141 | 0.078914 | 0.150634 | 0.953561 | 1.099773 | 9.072137 |
| Trial 5 | WO | UNADJ | 1.044529 | 0.025445 | 0.995831 | 1.095609 | 0.073703 | 1.000000 | 1.000000 | 0.000000 |
| Trial 5 | WO | LADJ | 1.033950 | 0.024017 | 0.987933 | 1.082111 | 0.150634 | 0.943903 | 1.122394 | 10.904724 |

## Elementary win/loss estimates

| Trial | Method | tau_w | tau_l | tau_tie |
| --- | --- | --- | --- | --- |
| Trial 5 | UNADJ | 0.215261 | 0.193481 | 0.591258 |
| Trial 5 | LADJ | 0.212542 | 0.195851 | 0.591607 |

## LADJ diagnostics

| Trial | Number_features | Gram_rank | Smallest_eigenvalue | Largest_eigenvalue | Condition_number | Scaled_condition_number | Max_centering_identity_error | Max_tau_formula_error | Max_scaled_fitted_value_error | Max_projection_normal_equation_error | Min_adjusted_covariance_eigenvalue |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Trial 5 | 10.000000 | 10.000000 | 85.364795 | 304355.654453 | 3565.353322 | 5.644213 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.173799 |

## LADJ point-estimate correction diagnostic

The feature-level LWR contributions below are a first-order decomposition of the exact elementary-probability correction. They are diagnostic and do not redefine the prespecified primary model.

| Feature | Raw_covariate | Standardized_mean_difference | Delta_tau_w | Delta_tau_l | First_order_delta_LWR |
| --- | --- | --- | --- | --- | --- |
| X_icu_0m1 | X_icu_0m | -0.036613 | -0.002099 | 0.002015 | -0.020167 |
| X_hospital_0m2 | X_hospital_0m | 0.028728 | -0.000585 | 0.000306 | -0.004301 |
| X_creat_0m | X_creat_0m | 0.019666 | 0.000390 | -0.000319 | 0.003458 |
| X_diabetes_elx_pmhx_0m1 | X_diabetes_elx_pmhx_0m | 0.030624 | -0.000320 | 0.000139 | -0.002204 |
| X_chf_elx_pmhx_0m1 | X_chf_elx_pmhx_0m | 0.035402 | -0.000123 | 0.000161 | -0.001404 |
| X_hospital_0m4 | X_hospital_0m | 0.008845 | 0.000021 | -0.000045 | 0.000329 |
| X_age_0m | X_age_0m | 0.011449 | -0.000030 | 0.000012 | -0.000201 |
| X_elx_score_pmhx_0m | X_elx_score_pmhx_0m | 0.010822 | 0.000022 | 0.000058 | -0.000195 |
| X_sex_0m1 | X_sex_0m | 0.005583 | -0.000032 | 0.000003 | -0.000167 |
| X_hospital_0m3 | X_hospital_0m | 0.017782 | 0.000039 | 0.000039 | -0.000025 |

## Leave-one-covariate-out diagnostic

These post-hoc fits diagnose the source of the adjustment and are not outcome-driven alternatives to the primary LADJ analysis.

| Model | Omitted_raw_covariate | LWR | ASE | CI_lower | CI_upper | P_value | Percent_variance_reduction | tau_w | tau_l | tau_tie |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Full primary LADJ |  | 0.081788 | 0.056910 | -0.029752 | 0.193329 | 0.150671 | 9.041501 | 0.212542 | 0.195851 | 0.591607 |
| Omit X_hospital_0m | X_hospital_0m | 0.084811 | 0.057068 | -0.027041 | 0.196663 | 0.137246 | 8.533024 | 0.212932 | 0.195617 | 0.591451 |
| Omit X_age_0m | X_age_0m | 0.081702 | 0.056912 | -0.029845 | 0.193248 | 0.151125 | 9.032317 | 0.212529 | 0.195856 | 0.591615 |
| Omit X_sex_0m | X_sex_0m | 0.081944 | 0.056912 | -0.029601 | 0.193490 | 0.149912 | 9.033007 | 0.212572 | 0.195848 | 0.591580 |
| Omit X_icu_0m | X_icu_0m | 0.100119 | 0.058946 | -0.015413 | 0.215651 | 0.089416 | 2.415289 | 0.214450 | 0.194019 | 0.591530 |
| Omit X_creat_0m | X_creat_0m | 0.079153 | 0.057108 | -0.032776 | 0.191083 | 0.165739 | 8.406326 | 0.212246 | 0.196094 | 0.591661 |
| Omit X_elx_score_pmhx_0m | X_elx_score_pmhx_0m | 0.081452 | 0.056890 | -0.030051 | 0.192955 | 0.152218 | 9.102750 | 0.212581 | 0.195953 | 0.591466 |
| Omit X_diabetes_elx_pmhx_0m | X_diabetes_elx_pmhx_0m | 0.083519 | 0.056920 | -0.028043 | 0.195081 | 0.142295 | 9.007100 | 0.212792 | 0.195742 | 0.591465 |
| Omit X_chf_elx_pmhx_0m | X_chf_elx_pmhx_0m | 0.082943 | 0.056924 | -0.028626 | 0.194512 | 0.145094 | 8.995519 | 0.212644 | 0.195718 | 0.591638 |

## Allocation-probability diagnostic

| Probability_definition | e | LWR | ASE |
| --- | --- | --- | --- |
| Known randomization probability | 0.500000 | 0.081788 | 0.056910 |
| Realized treatment proportion (diagnostic only) | 0.500395 | 0.081787 | 0.056909 |

## Quality control

| Check | Passed | Details |
| --- | --- | --- |
| All randomized subjects retained | TRUE | N = 5060 |
| Treatment has no missing values and both arms are present | TRUE | n1 = 2532 ; n0 = 2528 |
| All three 14-day outcomes are complete | TRUE | YS_death_14d 0; YS_dialysis_14d 0; YS_aki_progression_14d 0 |
| Outcome indicators are binary 0/1 | TRUE | Dictionary coding and observed values verified |
| Component union reproduces trial-defined 14-day composite | TRUE | YP_composite_14d |
| Every selected baseline covariate is complete | TRUE | X_hospital_0m, X_age_0m, X_sex_0m, X_icu_0m, X_creat_0m, X_elx_score_pmhx_0m, X_diabetes_elx_pmhx_0m, X_chf_elx_pmhx_0m |
| No post-randomization variable is used for LADJ | TRUE | All selected fields have dictionary role Baseline covariate |
| Pairwise win/loss/tie partition | TRUE | n1*n0 = 6400896 |
| LADJ Gram matrix full rank | TRUE | rank = 10 of 10 |
| LADJ projection normal equations | TRUE | maximum normal-equation residual = 1.073319e-11 |
| Adjusted elementary covariance is positive semidefinite | TRUE | minimum covariance eigenvalue = 1.737988e-01 |
| Exact transformations from elementary tau | TRUE | maximum error = 9.714451e-17 |
