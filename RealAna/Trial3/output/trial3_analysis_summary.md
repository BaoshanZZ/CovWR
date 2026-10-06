# Trial 3: 14-day hierarchical win analysis

## Verified mappings

| Role | Exact_variable | Coding_or_use |
| --- | --- | --- |
| Treatment | Treatment | A=1: Alert; A=0: Usual Care |
| 14-day death | YS_death_event | 1=death occurrence; 0=no death |
| 14-day dialysis | YS_dialysis_event | 1=dialysis occurrence; 0=no dialysis |
| 14-day AKI progression | YS_aki_progression_event | 1=AKI progression occurrence; 0=no progression |
| Hospital / study site | X_hospital_0m | Baseline LADJ covariate |
| Age category (fallback) | X_age_over_90_0m | Baseline LADJ covariate |
| Sex | X_sex_0m | Baseline LADJ covariate |
| ICU status | X_icu_0m | Baseline LADJ covariate |
| Elixhauser comorbidity score | X_elx_score_0m | Baseline LADJ covariate |
| Diabetes | X_diabetes_0m | Baseline LADJ covariate |
| Congestive heart failure | X_chf_0m | Baseline LADJ covariate |

## Analysis population

| Trial | Population | Role | N | Treatment_code_A | Known_randomization_probability |
| --- | --- | --- | --- | --- | --- |
| Trial 3 | All randomized | Total | 6030.0000 |  | 0.5000 |
| Trial 3 | Usual Care | Control | 2971.0000 | 0.0000 | 0.5000 |
| Trial 3 | Alert | Intervention | 3059.0000 | 1.0000 | 0.5000 |

## Outcome missingness

| Outcome | Variable | Missing |
| --- | --- | --- |
| 14-day death | YS_death_event | 0.0000 |
| 14-day dialysis | YS_dialysis_event | 0.0000 |
| 14-day AKI progression | YS_aki_progression_event | 0.0000 |

## Table A: outcomes by randomized arm

| Arm | Role | N | Death_n | Death_percent | Dialysis_n | Dialysis_percent | AKI_progression_n | AKI_progression_percent | Severity_death_n | Severity_death_percent | Severity_alive_dialysis_n | Severity_alive_dialysis_percent | Severity_alive_no_dialysis_AKI_n | Severity_alive_no_dialysis_AKI_percent | Severity_none_n | Severity_none_percent |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Usual Care | Control | 2971.0000 | 265.0000 | 8.9196 | 93.0000 | 3.1303 | 461.0000 | 15.5167 | 265.0000 | 8.9196 | 65.0000 | 2.1878 | 292.0000 | 9.8283 | 2349.0000 | 79.0643 |
| Alert | Intervention | 3059.0000 | 272.0000 | 8.8918 | 106.0000 | 3.4652 | 487.0000 | 15.9202 | 272.0000 | 8.8918 | 60.0000 | 1.9614 | 321.0000 | 10.4936 | 2406.0000 | 78.6532 |

## Unadjusted pairwise comparison

| Trial | Horizon_days | N | n1 | n0 | N_pairs | Treatment_wins | Treatment_losses | Final_ties | tau_w | tau_l | tau_tie | Decided_at_death_n | Decided_at_death_percent | Decided_at_dialysis_n | Decided_at_dialysis_percent | Decided_at_AKI_progression_n | Decided_at_AKI_progression_percent | Final_ties_percent |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Trial 3 | 14.0000 | 6030.0000 | 3059.0000 | 2971.0000 | 9088289.0000 | 1637306.0000 | 1669083.0000 | 5781900.0000 | 0.1802 | 0.1837 | 0.6362 | 1474587.0000 | 16.2251 | 352945.0000 | 3.8835 | 1478857.0000 | 16.2721 | 63.6192 |

## Table B: prespecified baseline covariate audit

| Clinical_covariate | Variable | Missing | Selected | Reason |
| --- | --- | --- | --- | --- |
| Hospital / study site | X_hospital_0m | 0.0000 | TRUE | Included: verified baseline and fully observed. |
| Age (continuous) | X_age_0m | 507.0000 | FALSE | Excluded: 507 missing value(s); no imputation or subject deletion. |
| Age category (fallback) | X_age_over_90_0m | 0.0000 | TRUE | Included as the fully observed prespecified age-category fallback. |
| Sex | X_sex_0m | 0.0000 | TRUE | Included: verified baseline and fully observed. |
| ICU status | X_icu_0m | 0.0000 | TRUE | Included: verified baseline and fully observed. |
| Baseline serum creatinine | X_creat_0m | 2.0000 | FALSE | Excluded: 2 missing value(s); no imputation or subject deletion. |
| Baseline AKI stage |  |  | FALSE | No true baseline AKI-stage variable is present in the dictionary. |
| Elixhauser comorbidity score | X_elx_score_0m | 0.0000 | TRUE | Included: verified baseline and fully observed. |
| Diabetes | X_diabetes_0m | 0.0000 | TRUE | Included: verified baseline and fully observed. |
| Congestive heart failure | X_chf_0m | 0.0000 | TRUE | Included: verified baseline and fully observed. |

**Final LADJ raw covariates:** X_hospital_0m, X_age_over_90_0m, X_sex_0m, X_icu_0m, X_elx_score_0m, X_diabetes_0m, X_chf_0m.

**Final LADJ design columns after reference dummy coding:** X_hospital_0m2, X_hospital_0m3, X_hospital_0m4, X_hospital_0m5, X_hospital_0m6, X_age_over_90_0m1, X_sex_0m1, X_icu_0m1, X_elx_score_0m, X_diabetes_0m1, X_chf_0m1.

## Table C: unadjusted and LADJ estimates

| Trial | Estimand | Method | Estimate | ASE | CI_lower | CI_upper | P_value | ASE_ratio_LADJ_to_UNADJ | Relative_efficiency_vs_UNADJ | Percent_variance_reduction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Trial 3 | NB | UNADJ | -0.003496 | 0.010602 | -0.024277 | 0.017284 | 0.741566 | 1.000000 | 1.000000 | 0.000000 |
| Trial 3 | NB | LADJ | -0.007910 | 0.010105 | -0.027716 | 0.011895 | 0.433738 | 0.953091 | 1.100857 | 9.161707 |
| Trial 3 | LWR | UNADJ | -0.019222 | 0.058284 | -0.133456 | 0.095012 | 0.741549 | 1.000000 | 1.000000 | 0.000000 |
| Trial 3 | LWR | LADJ | -0.043435 | 0.055492 | -0.152198 | 0.065328 | 0.433792 | 0.952102 | 1.103146 | 9.350192 |
| Trial 3 | WR | UNADJ | 0.980961 | 0.057174 | 0.875066 | 1.099672 | 0.741549 | 1.000000 | 1.000000 | 0.000000 |
| Trial 3 | WR | LADJ | 0.957495 | 0.053133 | 0.858819 | 1.067509 | 0.433792 | 0.929326 | 1.157881 | 13.635346 |
| Trial 3 | LWO | UNADJ | -0.006993 | 0.021205 | -0.048554 | 0.034568 | 0.741568 | 1.000000 | 1.000000 | 0.000000 |
| Trial 3 | LWO | LADJ | -0.015821 | 0.020211 | -0.055435 | 0.023793 | 0.433757 | 0.953139 | 1.100746 | 9.152559 |
| Trial 3 | WO | UNADJ | 0.993031 | 0.021057 | 0.952606 | 1.035173 | 0.741568 | 1.000000 | 1.000000 | 0.000000 |
| Trial 3 | WO | LADJ | 0.984303 | 0.019894 | 0.946074 | 1.024078 | 0.433757 | 0.944762 | 1.120354 | 10.742507 |

## Elementary win/loss estimates

| Trial | Method | tau_w | tau_l | tau_tie |
| --- | --- | --- | --- | --- |
| Trial 3 | UNADJ | 0.180156 | 0.183652 | 0.636192 |
| Trial 3 | LADJ | 0.178194 | 0.186105 | 0.635701 |

## LADJ diagnostics

| Trial | Number_features | Gram_rank | Smallest_eigenvalue | Largest_eigenvalue | Condition_number | Scaled_condition_number | Max_centering_identity_error | Max_tau_formula_error | Max_scaled_fitted_value_error | Max_projection_normal_equation_error | Min_adjusted_covariance_eigenvalue |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Trial 3 | 11.000000 | 11.000000 | 45.415167 | 12327.320711 | 271.436208 | 4.200346 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.173186 |

## LADJ point-estimate correction diagnostic

The feature-level LWR contributions below are a first-order decomposition of the exact elementary-probability correction. They are diagnostic and do not redefine the prespecified primary model.

| Feature | Raw_covariate | Standardized_mean_difference | Delta_tau_w | Delta_tau_l | First_order_delta_LWR |
| --- | --- | --- | --- | --- | --- |
| X_diabetes_0m1 | X_diabetes_0m | 0.036228 | -0.000708 | 0.000812 | -0.008350 |
| X_icu_0m1 | X_icu_0m | -0.019262 | -0.000701 | 0.000792 | -0.008205 |
| X_elx_score_0m | X_elx_score_0m | -0.013161 | -0.000634 | 0.000648 | -0.007047 |
| X_hospital_0m2 | X_hospital_0m | 0.026027 | -0.000310 | 0.000227 | -0.002954 |
| X_hospital_0m5 | X_hospital_0m | -0.011581 | 0.000152 | -0.000058 | 0.001163 |
| X_hospital_0m6 | X_hospital_0m | -0.014584 | 0.000142 | -0.000055 | 0.001089 |
| X_age_over_90_0m1 | X_age_over_90_0m | 0.016257 | 0.000053 | -0.000106 | 0.000876 |
| X_hospital_0m3 | X_hospital_0m | 0.008573 | -0.000043 | 0.000084 | -0.000695 |
| X_chf_0m1 | X_chf_0m | 0.003481 | -0.000054 | 0.000070 | -0.000682 |
| X_sex_0m1 | X_sex_0m | -0.026603 | 0.000111 | 0.000019 | 0.000513 |
| X_hospital_0m4 | X_hospital_0m | -0.028726 | 0.000030 | 0.000021 | 0.000052 |

## Leave-one-covariate-out diagnostic

These post-hoc fits diagnose the source of the adjustment and are not outcome-driven alternatives to the primary LADJ analysis.

| Model | Omitted_raw_covariate | LWR | ASE | CI_lower | CI_upper | P_value | Percent_variance_reduction | tau_w | tau_l | tau_tie |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Full primary LADJ |  | -0.043435 | 0.055492 | -0.152198 | 0.065328 | 0.433792 | 9.350192 | 0.178194 | 0.186105 | 0.635701 |
| Omit X_hospital_0m | X_hospital_0m | -0.042508 | 0.055645 | -0.151571 | 0.066554 | 0.444916 | 8.849452 | 0.178189 | 0.185927 | 0.635884 |
| Omit X_age_over_90_0m | X_age_over_90_0m | -0.044660 | 0.055498 | -0.153435 | 0.064114 | 0.420986 | 9.330351 | 0.178119 | 0.186254 | 0.635626 |
| Omit X_sex_0m | X_sex_0m | -0.044005 | 0.055517 | -0.152815 | 0.064806 | 0.427988 | 9.270083 | 0.178073 | 0.186084 | 0.635843 |
| Omit X_icu_0m | X_icu_0m | -0.037680 | 0.056460 | -0.148339 | 0.072980 | 0.504537 | 6.160424 | 0.178687 | 0.185548 | 0.635765 |
| Omit X_elx_score_0m | X_elx_score_0m | -0.031815 | 0.056572 | -0.142695 | 0.079065 | 0.573859 | 5.787014 | 0.179241 | 0.185035 | 0.635724 |
| Omit X_diabetes_0m | X_diabetes_0m | -0.034425 | 0.055774 | -0.143740 | 0.074891 | 0.537095 | 8.426350 | 0.178959 | 0.185227 | 0.635814 |
| Omit X_chf_0m | X_chf_0m | -0.041976 | 0.055658 | -0.151063 | 0.067112 | 0.450745 | 8.808049 | 0.178311 | 0.185955 | 0.635733 |

## Allocation-probability diagnostic

| Probability_definition | e | LWR | ASE |
| --- | --- | --- | --- |
| Known randomization probability | 0.500000 | -0.043435 | 0.055492 |
| Realized treatment proportion (diagnostic only) | 0.507297 | -0.043381 | 0.055539 |

## Quality control

| Check | Passed | Details |
| --- | --- | --- |
| All randomized subjects retained | TRUE | N = 6030 |
| Treatment has no missing values and both arms are present | TRUE | n1 = 3059 ; n0 = 2971 |
| All three 14-day outcomes are complete | TRUE | YS_death_event 0; YS_dialysis_event 0; YS_aki_progression_event 0 |
| Outcome indicators are binary 0/1 | TRUE | Dictionary coding and observed values verified |
| Component union reproduces trial-defined 14-day composite | TRUE | YP_composite_event |
| Every selected baseline covariate is complete | TRUE | X_hospital_0m, X_age_over_90_0m, X_sex_0m, X_icu_0m, X_elx_score_0m, X_diabetes_0m, X_chf_0m |
| No post-randomization variable is used for LADJ | TRUE | All selected fields have dictionary role Baseline covariate |
| Pairwise win/loss/tie partition | TRUE | n1*n0 = 9088289 |
| LADJ Gram matrix full rank | TRUE | rank = 11 of 11 |
| LADJ projection normal equations | TRUE | maximum normal-equation residual = 2.04729e-11 |
| Adjusted elementary covariance is positive semidefinite | TRUE | minimum covariance eigenvalue = 1.731858e-01 |
| Exact transformations from elementary tau | TRUE | maximum error = 2.220446e-16 |
