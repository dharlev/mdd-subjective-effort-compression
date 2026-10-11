# Analysis map

This file links each reported analysis to one public script and one explicit
input/output path. Supplementary sensitivities are kept beside the primary
analysis they qualify.

| Manuscript component | Public implementation | Model or method | Output |
|---|---|---|---|
| Table 1 and participant flow | `R/01_table1.R` | Welch/equal-variance t tests as labelled; Pearson chi-square; Cohen's d/phi | `outputs/tables/table1.csv` |
| Primary EP ratings | `R/02_primary_ep.R` | Hierarchical Gaussian rating model with correlated participant intercepts and effort slopes | `models/ep_primary.stan` |
| Age, sex, and success adjustment | `R/03_ep_sensitivities.R` | Joint hierarchical ANCOVA | `models/ep_covariates.stan` |
| Near-ceiling sensitivity | `R/03_ep_sensitivities.R` | Primary EP model after removing ratings above 98 | `outputs/ep_sensitivity/near_ceiling_summary.csv` |
| Actual execution sensitivities | `R/03_ep_sensitivities.R` | Rating model adjusted for force AUC, peak force, or mean force | `models/ep_flexible_covariates.stan` |
| Ordinal discrimination | `R/04_ordinal_piecewise.R` | Within-person ordered-pair score: 1 correct, 0 reversed, 0.5 tie; covariate-adjusted participant-level regression | `outputs/ordinal_piecewise/ordinal_summary.csv` |
| Piecewise shape | `R/04_ordinal_piecewise.R` | ML mixed models comparing one, two, or three adjacent-interval slopes by BIC | `outputs/ordinal_piecewise/piecewise_model_comparison.csv` |
| Dynamic observer | `R/05_dynamic_observer.R` | Dynamic Bayesian observer with group-specific prior and sensory uncertainty and shared learning rate | `models/dynamic_observer.stan` |
| Early force estimate | `R/06_effort_prior.R` | Joint hierarchical rating model with group-specific required-effort, early-estimate, and interaction coefficients | `models/force_prior.stan` |
| Execution-adjusted early estimate | `R/06_effort_prior.R` | Early-estimate model plus force AUC, success, age, and sex | `models/force_prior_execution_adjusted.stan` |
| Decision making | `R/07_decision_making.R` | Hierarchical Bernoulli-logit subjective-value model after catch-only exclusion | `models/decision_making.stan` |
| Symptoms and cross-task coupling | `R/08_clinical_cross_task.R` | Spearman association calculated for each posterior draw | `outputs/associations/` |
| Main/Supplementary empirical figures | `R/09_figures_tables.R` | Plots of exported posterior and descriptive summaries | `outputs/figures/` |

## Samples

- Recruited/contributor sample for Table 1: 94 controls, 50 participants with MDD.
- Primary EP analytic sample: 93 controls, 49 participants with MDD.
- Final decision-making model after the prespecified catch-trial exclusion:
  93 controls, 41 participants with MDD.

The code derives and checks these counts from input flags; it does not identify
or exclude participants by name.

## Standardisation

For EP models, rating, required effort, trial order, and continuous covariates
are standardised over the analytic observations/population specified in each
script. The original means and SDs are exported so coefficients can be returned
to the manuscript scale. Binary sex and success variables retain their documented
coding. The observer uses effort and rating on the 0-100 task scale.

For the decision-making model, reward is divided by the largest reward appearing
in the analytic offers. Effort is divided by 100 and squared in the subjective
value function. `choice_option2 = 1` means option 2 was selected; it does not mean
that the higher-effort option was necessarily selected.

