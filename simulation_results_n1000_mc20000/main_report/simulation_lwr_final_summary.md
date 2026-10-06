# Final LWR simulation summary: N = 1000 reference

## Design

- Reference sample size: N = 1000.
- Smaller-sample comparator: N = 500.
- Prognostic multipliers: 0.35 (weak), 0.70 (moderate), and 1.10 (strong).
- Alternative treatment-effect scale: 0.70.
- Methods: UNADJ and LADJ.
- Known randomization probability, centered no-intercept LADJ, closed-form Gram solution,
  and no clipping.
- Monte Carlo replications: 20,000 per scenario.
- Alternative truth evaluation: 10,000,000 independent treatment-control pairs.
- Raw replication-level results were not saved.

## Null performance

LWR type I error ranged from 4.90% to 5.11% for UNADJ and from 4.90% to 5.22% for LADJ.
The ASE/ESD ratio ranged from 0.993 to 1.003 for UNADJ and from 0.989 to 0.998 for LADJ.
Coverage remained close to 95% across all configurations.

## Alternative performance

The maximum absolute LWR bias was 0.00157 for UNADJ and 0.00195 for LADJ.

| Configuration | UNADJ power | LADJ power | LADJ PVR |
|---|---:|---:|---:|
| Reference | 43.24% | 53.35% | 21.50% |
| N = 500 | 24.37% | 30.18% | 21.73% |
| Allocation 1:2 | 39.37% | 48.43% | 22.66% |
| Weak prognosis | 55.29% | 58.74% | 7.76% |
| Strong prognosis | 30.60% | 45.62% | 37.54% |
| Effect heterogeneity | 43.60% | 52.53% | 20.74% |
| Endpoint-specific prognosis | 54.70% | 59.27% | 9.03% |
| Add noise covariates | 43.24% | 53.40% | 21.33% |
| Omit key covariate | 43.24% | 50.95% | 16.96% |

The results avoid power saturation. Strong prognosis produced the largest variance reduction,
whereas weak prognosis produced the smallest gain. Adding irrelevant covariates had little
effect relative to the reference. Omitting X3 reduced both PVR and power.
