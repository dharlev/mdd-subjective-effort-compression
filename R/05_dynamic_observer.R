prepare_dynamic_observer <- function(inputs) {
  participants <- inputs$participants %>%
    filter(ep_primary_include == 1) %>%
    arrange(diagnosis, participant_id) %>%
    mutate(subject_index = row_number(), group_id = if_else(diagnosis == "control", 1L, 2L))
  trials <- inputs$ep %>%
    inner_join(participants %>% select(participant_id, subject_index),
               by = "participant_id") %>%
    filter(include_primary == 1, is.finite(required_effort_pct),
           is.finite(rating), is.finite(trial)) %>%
    group_by(subject_index) %>%
    arrange(trial, .by_group = TRUE) %>%
    mutate(trial_z_within = z_safe(trial)) %>%
    ungroup() %>%
    arrange(subject_index, trial)
  # Explicit row bounds avoid relying on trial labels or equal trial counts.
  run_lengths <- rle(trials$subject_index)$lengths
  subj_end <- cumsum(run_lengths)
  subj_start <- c(1L, head(subj_end, -1L) + 1L)

  stan_data <- list(
    N = nrow(trials), N_subj = nrow(participants), G = 2L,
    group = participants$group_id,
    subj_start = as.integer(subj_start), subj_end = as.integer(subj_end),
    effort_pct = trials$required_effort_pct,
    trial_z = trials$trial_z_within,
    y = trials$rating,
    beta_trial_fixed = 0,
    sigma_response_fixed = 0.001
  )
  list(participants = participants, trials = trials, stan_data = stan_data)
}

summarise_observer <- function(fit, prepared, out) {
  prior <- draw_matrix(fit, "sigma_prior_group")
  sensory <- draw_matrix(fit, "sigma_sensory_group")
  weight <- draw_matrix(fit, "w_prior_group")
  initial <- draw_matrix(fit, "mu_prior0_group")
  lambda <- as.numeric(draw_matrix(fit, "lambda"))

  quantities <- list(
    control_prior_uncertainty = prior[, 1], MDD_prior_uncertainty = prior[, 2],
    MDD_minus_control_prior_uncertainty = prior[, 2] - prior[, 1],
    control_sensory_uncertainty = sensory[, 1], MDD_sensory_uncertainty = sensory[, 2],
    MDD_minus_control_sensory_uncertainty = sensory[, 2] - sensory[, 1],
    control_effective_prior_weight = weight[, 1], MDD_effective_prior_weight = weight[, 2],
    MDD_minus_control_effective_prior_weight = weight[, 2] - weight[, 1],
    control_initial_prior_mean = initial[, 1], MDD_initial_prior_mean = initial[, 2],
    MDD_minus_control_initial_prior_mean = initial[, 2] - initial[, 1],
    shared_learning_rate = lambda
  )
  summary <- imap_dfr(quantities, ~ summarise_draws_vector(.x, .y, "task scale"))
  write_csv(summary, file.path(out, "observer_summary.csv"))

  y_rep <- draw_matrix(fit, "y_rep")
  group <- prepared$participants$diagnosis[prepared$trials$subject_index]
  ppc <- expand_grid(
    diagnosis = DIAGNOSIS_LEVELS,
    required_effort_pct = sort(unique(prepared$trials$required_effort_pct))
  ) %>%
    rowwise() %>%
    mutate(
      draw_values = list(rowMeans(y_rep[, group == diagnosis &
                                          prepared$trials$required_effort_pct == required_effort_pct,
                                        drop = FALSE])),
      observed = mean(prepared$trials$rating[group == diagnosis &
                                               prepared$trials$required_effort_pct == required_effort_pct]),
      predicted = mean(draw_values),
      q025 = quantile(draw_values, 0.025),
      q975 = quantile(draw_values, 0.975)
    ) %>%
    ungroup() %>%
    select(-draw_values)
  write_csv(ppc, file.path(out, "observer_posterior_predictive.csv"))
  list(summary = summary, ppc = ppc)
}

run_dynamic_observer <- function(inputs) {
  out <- file.path(OUTPUT_DIR, "dynamic_observer")
  prepared <- prepare_dynamic_observer(inputs)
  fit <- sample_stan("dynamic_observer.stan", prepared$stan_data, out,
                     "dynamic_observer", seed = 20260617L)
  summary <- summarise_observer(fit, prepared, out)
  list(fit = fit, prepared = prepared, summary = summary)
}
