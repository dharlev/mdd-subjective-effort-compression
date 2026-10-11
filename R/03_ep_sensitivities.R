prepare_ep_covariate_model <- function(inputs, execution_metric = NULL) {
  participants <- inputs$participants %>%
    filter(ep_primary_include == 1) %>%
    arrange(source_group, participant_id) %>%
    mutate(
      subject_index = row_number(),
      group_id = as.integer(source_group),
      age_z = z_safe(age_years),
      sex_male = as.numeric(sex == "male")
    )
  trials <- inputs$ep %>%
    inner_join(participants %>%
                 select(participant_id, subject_index, age_z, sex_male),
               by = "participant_id") %>%
    filter(include_primary == 1, is_practice == 0) %>%
    arrange(subject_index, trial) %>%
    mutate(
      rating_z = z_safe(rating),
      effort_z = z_safe(required_effort_pct),
      trial_z = z_safe(trial),
      success_z = z_safe(success)
    )

  metric_z <- NULL
  if (!is.null(execution_metric)) {
    if (!execution_metric %in% c("force_auc", "peak_force", "mean_force")) {
      stop("Unknown execution metric.")
    }
    if (any(!is.finite(trials[[execution_metric]]))) {
      stop("Missing values in execution metric: ", execution_metric)
    }
    metric_z <- z_safe(trials[[execution_metric]])
  }

  X_subj <- as.matrix(participants[, c("age_z", "sex_male")])
  X_obs <- cbind(trial_z = trials$trial_z, success_z = trials$success_z)
  if (!is.null(metric_z)) X_obs <- cbind(X_obs, execution_metric_z = metric_z)

  if (is.null(execution_metric)) {
    stan_data <- list(
      N_subj = nrow(participants), N_obs = nrow(trials),
      N_group = length(GROUP_LEVELS), K_subj = ncol(X_subj), K_obs = ncol(X_obs),
      subj = trials$subject_index, subj_group = participants$group_id,
      rating_z = trials$rating_z, effort_z = trials$effort_z,
      X_subj = X_subj, X_obs = X_obs
    )
    model_file <- "ep_covariates.stan"
  } else {
    stan_data <- list(
      N_subj = nrow(participants), N_obs = nrow(trials),
      N_group = length(GROUP_LEVELS), K_alpha = ncol(X_subj),
      K_beta = ncol(X_subj), K_obs = ncol(X_obs),
      subj = trials$subject_index, subj_group = participants$group_id,
      rating_z = trials$rating_z, effort_z = trials$effort_z,
      X_alpha = X_subj, X_beta = X_subj, X_obs = X_obs
    )
    model_file <- "ep_flexible_covariates.stan"
  }
  list(participants = participants, trials = trials, stan_data = stan_data,
       model_file = model_file, execution_metric = execution_metric)
}

summarise_ep_covariate_fit <- function(fit, prepared, output_file) {
  participants <- prepared$participants
  weights <- control_weights(participants, "ep_primary_include")
  mu_alpha <- draw_matrix(fit, "mu_alpha_group")
  mu_beta <- draw_matrix(fit, "mu_beta_group")
  alpha <- draw_matrix(fit, "alpha")
  beta <- draw_matrix(fit, "beta")
  idx_c <- which(participants$diagnosis == "control")
  idx_m <- which(participants$diagnosis == "MDD")

  fixed <- list(
    control_intercept = weighted_control_draw(mu_alpha, weights),
    MDD_intercept = mu_alpha[, 3],
    control_discriminability = weighted_control_draw(mu_beta, weights),
    MDD_discriminability = mu_beta[, 3]
  )
  fixed$MDD_minus_control_intercept <- fixed$MDD_intercept - fixed$control_intercept
  fixed$MDD_minus_control_discriminability <-
    fixed$MDD_discriminability - fixed$control_discriminability

  observed_covariate_average <- list(
    control_intercept = rowMeans(alpha[, idx_c, drop = FALSE]),
    MDD_intercept = rowMeans(alpha[, idx_m, drop = FALSE]),
    control_discriminability = rowMeans(beta[, idx_c, drop = FALSE]),
    MDD_discriminability = rowMeans(beta[, idx_m, drop = FALSE])
  )
  observed_covariate_average$MDD_minus_control_intercept <-
    observed_covariate_average$MDD_intercept - observed_covariate_average$control_intercept
  observed_covariate_average$MDD_minus_control_discriminability <-
    observed_covariate_average$MDD_discriminability -
    observed_covariate_average$control_discriminability

  output <- bind_rows(
    imap_dfr(fixed, ~ summarise_draws_vector(.x, .y) %>%
               mutate(estimand = "group_fixed_effect_at_mean_age_reference_sex", .before = 1)),
    imap_dfr(observed_covariate_average, ~ summarise_draws_vector(.x, .y) %>%
               mutate(estimand = "own_group_observed_covariate_average", .before = 1))
  )
  write_csv(output, output_file)
  output
}

run_ep_sensitivities <- function(inputs) {
  out <- file.path(OUTPUT_DIR, "ep_sensitivity")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  adjusted_data <- prepare_ep_covariate_model(inputs)
  adjusted_fit <- sample_stan(adjusted_data$model_file, adjusted_data$stan_data,
                              file.path(out, "age_sex_success"), "ep_age_sex_success",
                              seed = 20260726L)
  adjusted_summary <- summarise_ep_covariate_fit(
    adjusted_fit, adjusted_data, file.path(out, "age_sex_success_summary.csv")
  )

  ceiling_data <- prepare_primary_ep(inputs, ceiling_threshold = 98)
  ceiling_fit <- sample_stan("ep_primary.stan", ceiling_data$stan_data,
                            file.path(out, "near_ceiling"), "ep_near_ceiling",
                            seed = 20260723L)
  ceiling_summary <- summarise_primary_ep(
    ceiling_fit, ceiling_data, file.path(out, "near_ceiling"), "near_ceiling"
  )
  ceiling_removal <- inputs$ep %>%
    inner_join(inputs$participants %>%
                 filter(ep_primary_include == 1) %>%
                 select(participant_id, diagnosis), by = "participant_id") %>%
    filter(include_primary == 1, is_practice == 0) %>%
    group_by(diagnosis) %>%
    summarise(total = n(), removed = sum(rating > 98),
              percent_removed = 100 * removed / total, .groups = "drop")
  write_csv(ceiling_removal, file.path(out, "near_ceiling_trials_removed.csv"))

  execution <- map(c("force_auc", "peak_force", "mean_force"), function(metric) {
    prepared <- prepare_ep_covariate_model(inputs, metric)
    metric_out <- file.path(out, paste0("execution_", metric))
    fit <- sample_stan(prepared$model_file, prepared$stan_data, metric_out,
                       paste0("ep_execution_", metric), seed = 20260726L)
    summary <- summarise_ep_covariate_fit(
      fit, prepared, file.path(metric_out, "summary.csv")
    ) %>% mutate(execution_metric = metric, .before = 1)
    list(metric = metric, fit = fit, prepared = prepared, summary = summary)
  })
  execution_summary <- bind_rows(map(execution, "summary"))
  write_csv(execution_summary, file.path(out, "execution_metric_comparison.csv"))

  list(
    adjusted = list(fit = adjusted_fit, prepared = adjusted_data,
                    summary = adjusted_summary),
    ceiling = list(fit = ceiling_fit, prepared = ceiling_data,
                   summary = ceiling_summary, removal = ceiling_removal),
    execution = execution,
    execution_summary = execution_summary
  )
}
