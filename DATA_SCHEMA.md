# Input data contract

The public pipeline requires four deidentified CSV files in `data/`. Data are
not included in the repository. Column names and coding must match exactly.

## Privacy requirement

`participant_id` must be a release-specific anonymous identifier matching
`P####` (for example, `P0001`). Do not use study IDs, initials, dates, folder
names, recruitment-site labels, free text, or any direct/quasi identifier.

## `participants.csv`

One row per person who contributed to any analysis.

| Column | Coding |
|---|---|
| `participant_id` | Anonymous `P####` identifier |
| `diagnosis` | `control` or `MDD` |
| `control_cohort` | `cohort_1`, `cohort_2`, or blank for MDD |
| `age_years` | Numeric |
| `sex` | `female` or `male` |
| `phq9_total` | Numeric; blank when unavailable |
| `aes_total` | Numeric; blank when unavailable |
| `ami_total` | Numeric; blank when unavailable |
| `table1_include` | 0/1 |
| `ep_primary_include` | 0/1 |
| `dm_task_complete` | 0/1 |
| `dm_model_include` | 0/1; catch-only final sample |

`dm_model_include` must be derived only from the prespecified catch criterion:
exclude a participant when the non-dominant option was chosen on both catch
trials. Do not exclude extreme or near-deterministic ordinary choice patterns.

## `ep_trials.csv`

One row per effort-psychophysics trial.

| Column | Coding/unit |
|---|---|
| `participant_id` | Links to `participants.csv` |
| `trial` | Within-task trial order |
| `is_practice` | 0/1 |
| `required_effort_pct` | 10, 30, 50, or 70 (%MVC) |
| `rating` | Subjective effort rating, 0-100 |
| `success` | 0/1 |
| `force_auc` | Positive baseline-corrected force AUC, MVC-normalised |
| `peak_force` | Peak baseline-corrected force, MVC-normalised |
| `mean_force` | Mean baseline-corrected force, MVC-normalised |
| `early_force_mvc` | Extracted early force estimate, proportion MVC |
| `early_force_mvc_20` | Optional 20-sample-kernel sensitivity estimate |
| `early_force_time_samples` | Sample index of estimate, relative to force-trace onset |
| `include_primary` | 0/1 |
| `include_early_force` | 0/1 |

Force traces were sampled at 250 Hz. For the early estimate, force was
baseline-corrected using the first 10 samples; the derivative was smoothed with
a centred 40-sample (160 ms) kernel; the anchor was the maximum smoothed
derivative; the estimate was the first later sample with force above 0.15 MVC
and absolute derivative below 0.001. Timing outliers beyond three participant
SDs were excluded. A 20-sample (80 ms) extraction was a sensitivity analysis,
not the primary estimate.

For execution metrics, force was baseline-corrected using the first 10 samples
and MVC-normalised. `force_auc` is the trapezoidal integral of positive
baseline-corrected force across the full trace; `peak_force` is its maximum;
`mean_force` is AUC divided by trace duration.

## `dm_trials.csv`

One row per decision-making trial.

| Column | Coding/unit |
|---|---|
| `participant_id` | Links to `participants.csv` |
| `trial` | Within-task trial order |
| `is_practice` | 0/1 |
| `is_catch` | 0/1 |
| `choice_option2` | 1 if option 2 selected, otherwise 0 |
| `reward1`, `reward2` | Offered points |
| `effort1_pct`, `effort2_pct` | Offered effort as %MVC |
| `success` | 0/1, when applicable |

## `medications.csv` (optional)

One row per listed medication for participants with MDD.

| Column | Coding |
|---|---|
| `participant_id` | Links to `participants.csv` |
| `medication_class` | `antidepressant`, `antipsychotic`, `sedative`, or `other` |
| `medication_name` | Generic medication name only |
| `dose_value` | Numeric dose; blank if unavailable |
| `dose_unit` | Normalised unit such as `mg` |
| `prn` | 0/1 |

Relative burden is the dose divided by the median positive dose for the same
medication and unit within the MDD sample, summed overall and within class.
Entries without a numeric dose contribute to counts/class indicators but not to
continuous burden. PRN entries are included as listed without assumed frequency.
