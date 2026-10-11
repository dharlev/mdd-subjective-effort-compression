spearman_draws_vector <- function(parameter_draws, outcome) {
  apply(parameter_draws, 1, function(x) cor(x, outcome, method = "spearman",
                                           use = "complete.obs"))
}

spearman_draws_paired <- function(x_draws, y_draws) {
  stopifnot(nrow(x_draws) == nrow(y_draws), ncol(x_draws) == ncol(y_draws))
  vapply(seq_len(nrow(x_draws)), function(i) {
    cor(x_draws[i, ], y_draws[i, ], method = "spearman", use = "complete.obs")
  }, numeric(1))
}

prepare_symptom_trial_model <- function(inputs) {
  participants <- inputs$participants %>%
    filter(ep_primary_include == 1, diagnosis == "MDD", is.finite(phq9_total)) %>%
    arrange(participant_id) %>%
    mutate(subject_index = row_number(), phq_z = z_safe(phq9_total),
           age_z = z_safe(age_years), sex_male = as.numeric(sex == "male"))
  trials <- inputs$ep %>%
    inner_join(participants %>%
                 select(participant_id, subject_index, phq_z, age_z, sex_male),
               by = "participant_id") %>%
    filter(include_primary == 1, is_practice == 0) %>%
    arrange(subject_index, trial) %>%
    mutate(rating_z = z_safe(rating), effort_z = z_safe(required_effort_pct),
           trial_z = z_safe(trial))
  list(
    participants = participants, trials = trials,
    stan_data = list(
      N_subj = nrow(participants), N_obs = nrow(trials),
      subj = trials$subject_index, rating_z = trials$rating_z,
      effort_z = trials$effort_z, trial_z = trials$trial_z,
      phq_z = trials$phq_z, age_z = trials$age_z,
      sex_male = trials$sex_male
    )
  )
}

run_symptom_trial_model <- function(inputs, out) {
  prepared <- prepare_symptom_trial_model(inputs)
  fit <- sample_stan("ep_symptom.stan", prepared$stan_data, out,
                     "ep_MDD_PHQ_trial_model")
  variables <- c("beta_phq", "beta_phq_effort", "beta_age", "beta_sex", "beta_trial")
  summary <- map_dfr(variables, function(v) {
    summarise_draws_vector(as.numeric(draw_matrix(fit, v)), v)
  })
  write_csv(summary, file.path(out, "ep_MDD_PHQ_trial_summary.csv"))
  list(fit = fit, prepared = prepared, summary = summary)
}

construct_medication_burden <- function(medications) {
  if (is.null(medications)) return(NULL)
  med <- medications %>%
    mutate(
      medication_name = tolower(trimws(medication_name)),
      dose_unit = tolower(trimws(dose_unit)),
      positive_dose = if_else(is.finite(dose_value) & dose_value > 0,
                              dose_value, NA_real_)
    ) %>%
    group_by(medication_name, dose_unit) %>%
    mutate(reference_median = median(positive_dose, na.rm = TRUE),
           relative_dose = if_else(is.finite(positive_dose) &
                                     is.finite(reference_median) & reference_median > 0,
                                   positive_dose / reference_median, 0)) %>%
    ungroup()
  med %>%
    group_by(participant_id) %>%
    summarise(
      medication_count = n(),
      total_burden = sum(relative_dose),
      antidepressant_burden = sum(relative_dose[medication_class == "antidepressant"]),
      antipsychotic_burden = sum(relative_dose[medication_class == "antipsychotic"]),
      sedative_burden = sum(relative_dose[medication_class == "sedative"]),
      antidepressant = as.numeric(any(medication_class == "antidepressant")),
      antipsychotic = as.numeric(any(medication_class == "antipsychotic")),
      sedative = as.numeric(any(medication_class == "sedative")),
      .groups = "drop"
    )
}

medication_associations <- function(inputs, ep_primary, out) {
  burden <- construct_medication_burden(inputs$medications)
  if (is.null(burden)) return(NULL)
  p <- ep_primary$prepared$participants
  mdd_idx <- which(p$diagnosis == "MDD")
  mdd <- p[mdd_idx, ] %>%
    select(participant_id, age_years, sex) %>%
    left_join(burden, by = "participant_id") %>%
    mutate(across(c(medication_count, total_burden, antidepressant_burden,
                    antipsychotic_burden, sedative_burden,
                    antidepressant, antipsychotic, sedative),
                  ~ replace_na(.x, 0)))
  predictors <- c("medication_count", "total_burden", "antidepressant_burden",
                  "antipsychotic_burden", "sedative_burden")
  alpha <- ep_primary$summary$alpha_draws[, mdd_idx, drop = FALSE]
  beta <- ep_primary$summary$beta_draws[, mdd_idx, drop = FALSE]
  output <- map_dfr(predictors, function(predictor) {
    x <- z_safe(mdd[[predictor]])
    if (!all(is.finite(x))) return(NULL)
    slope_beta <- apply(beta, 1, function(y) coef(lm(y ~ x))[2])
    slope_alpha <- apply(alpha, 1, function(y) coef(lm(y ~ x))[2])
    bind_rows(
      summarise_draws_vector(slope_beta, "effort_discriminability"),
      summarise_draws_vector(slope_alpha, "rating_intercept")
    ) %>% mutate(medication_predictor = predictor, .before = 1)
  })
  write_csv(output, file.path(out, "medication_within_MDD_associations.csv"))
  write_csv(burden, file.path(out, "medication_burden_deidentified.csv"))
  output
}

run_clinical_cross_task <- function(inputs, ep_primary, dm) {
  out <- file.path(OUTPUT_DIR, "associations")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  ep_p <- ep_primary$prepared$participants
  alpha <- ep_primary$summary$alpha_draws
  beta <- ep_primary$summary$beta_draws

  clinical <- map_dfr(DIAGNOSIS_LEVELS, function(group) {
    idx <- which(ep_p$diagnosis == group)
    phq <- ep_p$phq9_total[idx]
    aes <- ep_p$aes_total[idx]
    bind_rows(
      summarise_draws_vector(spearman_draws_vector(beta[, idx, drop = FALSE], phq),
                             "discriminability_vs_PHQ", "Spearman rho"),
      summarise_draws_vector(spearman_draws_vector(alpha[, idx, drop = FALSE], phq),
                             "intercept_vs_PHQ", "Spearman rho"),
      summarise_draws_vector(spearman_draws_vector(beta[, idx, drop = FALSE], aes),
                             "discriminability_vs_AES", "Spearman rho"),
      summarise_draws_vector(spearman_draws_vector(alpha[, idx, drop = FALSE], aes),
                             "intercept_vs_AES", "Spearman rho")
    ) %>% mutate(diagnosis = group, .before = 1)
  })
  write_csv(clinical, file.path(out, "clinical_associations.csv"))

  dm_p <- dm$prepared$participants
  common <- intersect(ep_p$participant_id, dm_p$participant_id)
  ep_order <- match(common, ep_p$participant_id)
  dm_order <- match(common, dm_p$participant_id)
  common_group <- as.character(ep_p$diagnosis[ep_order])
  dm_effort <- dm$summary$effort_draws[, dm_order, drop = FALSE]
  cross_task <- map_dfr(DIAGNOSIS_LEVELS, function(group) {
    keep <- common_group == group
    slope_rho <- spearman_draws_paired(beta[, ep_order[keep], drop = FALSE],
                                      dm_effort[, keep, drop = FALSE])
    intercept_rho <- spearman_draws_paired(alpha[, ep_order[keep], drop = FALSE],
                                          dm_effort[, keep, drop = FALSE])
    bind_rows(
      summarise_draws_vector(slope_rho, "EP_discriminability_vs_DM_effort", "Spearman rho"),
      summarise_draws_vector(intercept_rho, "EP_intercept_vs_DM_effort", "Spearman rho"),
      summarise_draws_vector(slope_rho - intercept_rho,
                             "discriminability_minus_intercept_specificity", "rho difference")
    ) %>% mutate(diagnosis = group, .before = 1)
  })
  write_csv(cross_task, file.path(out, "cross_task_associations.csv"))

  trial_model <- run_symptom_trial_model(inputs, file.path(out, "PHQ_trial_model"))
  medication <- medication_associations(inputs, ep_primary, out)
  list(clinical = clinical, cross_task = cross_task,
       trial_model = trial_model, medication = medication)
}

