# Not just harder: subjective effort compression in major depressive disorder

Reproducible R and Stan code for the analyses, numerical tables, and empirical
figures reported in:

> Harlev D, Wolpe N, et al. *Not just harder: A compression of perceived effort
> in major depressive disorder is associated with symptom burden and altered
> decision-making.*

The repository is intentionally organised as a publication pipeline rather than
as a record of exploratory development. It starts from four deidentified,
analysis-ready tables, applies the manuscript exclusions, fits the reported
models, and writes auditable numerical summaries before drawing figures.

Participant-level data are not stored in this public repository. Access is
governed by the study's ethics approval and Data Availability statement. The
required columns, coding, units, and exclusion flags are specified in
[`DATA_SCHEMA.md`](DATA_SCHEMA.md).

## Quick start

1. Install R (the final analyses were checked with R 4.4.2), CmdStan, and the
   packages listed in `R/00_setup.R`. `Rscript install_packages.R` provides a
   convenience installer.
2. Place the deidentified CSV files described in `DATA_SCHEMA.md` in `data/`.
3. From the repository root, run:

```sh
Rscript run_all.R
```

The complete Bayesian pipeline uses four chains, 1,000 warm-up iterations and
4,000 post-warm-up iterations per chain. It can therefore take several hours.
For a code and data-contract check without fitting models, run:

```sh
Rscript validate_release.R
```

## Pipeline

| Stage | Script | Main products |
|---|---|---|
| 1 | `R/01_table1.R` | Participant flow and Table 1 |
| 2 | `R/02_primary_ep.R` | Primary effort-rating intercept and discriminability analyses |
| 3 | `R/03_ep_sensitivities.R` | Covariate, near-ceiling, and actual-execution sensitivities |
| 4 | `R/04_ordinal_piecewise.R` | Ordinal discrimination and piecewise shape analyses |
| 5 | `R/05_dynamic_observer.R` | Dynamic observer fit and posterior predictive checks |
| 6 | `R/06_effort_prior.R` | Early force estimate to rating coupling |
| 7 | `R/07_decision_making.R` | Catch-only decision-making model and recovery inputs |
| 8 | `R/08_clinical_cross_task.R` | Symptom and cross-task posterior-draw associations |
| 9 | `R/09_figures_tables.R` | Manuscript tables and empirical figures |

[`ANALYSIS_MAP.md`](ANALYSIS_MAP.md) maps each manuscript result to its source
script and model. [`EXPECTED_RESULTS.md`](EXPECTED_RESULTS.md) provides rounded
reference values for checking a reproduction.

## Important estimands

The primary effort-psychophysics fit reports two related quantities explicitly:

- the **finite-sample participant-average contrast**, obtained by averaging the
  fitted participant intercepts and slopes within the enrolled analytic sample;
- the **hierarchical population contrast**, obtained from the group-level mean
  parameters.

The manuscript's reported contrast of approximately `-0.072` is the first
estimand. Both are exported with unambiguous labels so they cannot be mistaken
for one another.

The two healthy recruitment cohorts are represented separately in all relevant
hierarchical models. A single control estimate is formed with their fixed EP or
DM analytic sample-size weights, as appropriate.

## Reproducibility and privacy

- Random seeds are fixed in `R/00_setup.R`.
- Every fitted model exports convergence diagnostics and posterior summaries.
- All paths are relative to the repository root.
- No names, source-study IDs, dates of participation, local paths, or raw force
  traces are included.
- The pipeline rejects participant identifiers that do not follow the anonymous
  `P####` format.
- Generated outputs, fitted model objects, and private input CSV files are
  excluded by `.gitignore`.

The conceptual task/hypothesis artwork is supplied in `assets/`. Empirical
panels are regenerated from model summaries by `R/09_figures_tables.R`.

## Software

The main dependencies are `cmdstanr`, `posterior`, `dplyr`, `tidyr`, `readr`,
`purrr`, `ggplot2`, `patchwork`, `nlme`, `MASS`, `loo`, and `scales`. Exact
installed versions are recorded in `outputs/software_versions.csv` and
`outputs/sessionInfo.txt` on every run.

## Licence

Code is released under the MIT License. The data are not covered by this code
licence.
