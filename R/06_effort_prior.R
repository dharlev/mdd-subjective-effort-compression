prepare_effort_prior <- function(inputs, prior_column = "early_force_mvc") {
  if (!prior_column %in% names(inputs$ep)) stop("Missing early-force column: ", prior_column)
  participants <- inputs$participants %>%
    filter(ep_primary_include == 1) %>%
    semi_join(inputs$ep %>% filter(include_early_force == 1,
                                   is.finite(.data[[prior_column]])),
              by = "participant_id") %>%
    arrange(diagnosis, participant_id) %>%
    mutate(subject_index = row_number(), group_id = if_else(diagnosis == "control", 1L, 2L))
  trials <- inputs$ep %>%
    inner_join(participants %>% select(participant_id, subject_index),
               by = "participant_id") %>%
    filter(include_early_force == 1, is_practice == 0,
           is.finite(.data[[prior_column]]), is.finite(rating)) %>%
    group_by(subject_index) %>%
    arrange(trial, .by_group = TRUE) %>%
    mutate(trial_z = z_safe(trial)) %>%
    ungroup() %>%
    mutate(
      req_force = required_effort_pct / 100,
      prior_force = .data[[prior_column]]
    ) %>%
    arrange(subject_index, trial)
  stan_data <- list(
    N = nrow(trials), N_subj = nrow(participants), G = 2L,
    subj = trials$subject_index, group = participants$group_id,
    req_force = trials$req_force, prior_force = trials$prior_force,
    trial_z = trials$trial_z, rating_obs = trials$rating
  )
  list(participants = participants, trials = trials, stan_data = stan_data,
       prior_column = prior_column)
}

summarise_effort_prior <- function(fit, prepared, out, prefix = "primary_40_sample") {
  effort <- draw_matrix(fit, "beta_effort_group")
  prior <- draw_matrix(fit, "beta_prior_group")
  interaction <- draw_matrix(fit, "beta_interaction_group")
  d <- prepared$trials
  mu_req <- mean(d$req_force); sd_req <- sd(d$req_force)
  mu_prior <- mean(d$prior_force); sd_prior <- sd(d$prior_force)
  sd_rating <- sd(d$rating)

  effort_z <- sweep(effort, 1, interaction[, 1] * 0, "+")
  effort_z[, 1] <- (effort[, 1] + interaction[, 1] * mu_prior) * sd_req / sd_rating
  effort_z[, 2] <- (effort[, 2] + interaction[, 2] * mu_prior) * sd_req / sd_rating
  prior_z <- prior
  prior_z[, 1] <- (prior[, 1] + interaction[, 1] * mu_req) * sd_prior / sd_rating
  prior_z[, 2] <- (prior[, 2] + interaction[, 2] * mu_req) * sd_prior / sd_rating
  interaction_z <- interaction * sd_req * sd_prior / sd_rating

  quantities <- list(
    control_required_effort = effort_z[, 1], MDD_required_effort = effort_z[, 2],
    MDD_minus_control_required_effort = effort_z[, 2] - effort_z[, 1],
    control_early_estimate = prior_z[, 1], MDD_early_estimate = prior_z[, 2],
    MDD_minus_control_early_estimate = prior_z[, 2] - prior_z[, 1],
    control_early_estimate_by_effort = interaction_z[, 1],
    MDD_early_estimate_by_effort = interaction_z[, 2],
    MDD_minus_control_early_estimate_by_effort = interaction_z[, 2] - interaction_z[, 1]
  )
  summary <- imap_dfr(quantities, ~ summarise_draws_vector(.x, .y)) %>%
    mutate(extraction = prefix, .before = 1)
  write_csv(summary, file.path(out, paste0(prefix, "_summary.csv")))

  levels <- c(0.10, 0.30, 0.50, 0.70)
  simple <- map_dfr(levels, function(x) {
    control <- (prior[, 1] + interaction[, 1] * x) * sd_prior / sd_rating
    mdd <- (prior[, 2] + interaction[, 2] * x) * sd_prior / sd_rating
    bind_rows(
      summarise_draws_vector(control, "control_early_estimate_coupling"),
      summarise_draws_vector(mdd, "MDD_early_estimate_coupling"),
      summarise_draws_vector(mdd - control, "MDD_minus_control_early_estimate_coupling")
    ) %>% mutate(required_effort_pct = x * 100, .before = 1)
  })
  write_csv(simple, file.path(out, paste0(prefix, "_simple_coupling.csv")))
  list(summary = summary, simple = simple, standardised = list(
    effort = effort_z, prior = prior_z, interaction = interaction_z
  ))
}

prepare_effort_prior_execution <- function(prepared) {
  p <- prepared$participants %>%
    left_join(
      prepared$trials %>% distinct(participant_id) %>% mutate(present = TRUE),
      by = "participant_id"
    )
  d <- prepared$trials %>%
    left_join(p %>% select(participant_id, age_years, sex), by = "participant_id") %>%
    mutate(
      req_force_z = z_safe(req_force), prior_force_z = z_safe(prior_force),
      rating_z = z_safe(rating), auc_force_z = z_safe(force_auc),
      success_z = z_safe(success), age_z = z_safe(age_years),
      sex_male = as.numeric(sex == "male")
    )
  list(
    participants = p, trials = d,
    stan_data = list(
      N = nrow(d), N_subj = nrow(p), G = 2L,
      subj = d$subject_index, group = p$group_id,
      req_force_z = d$req_force_z, prior_force_z = d$prior_force_z,
      trial_z = d$trial_z, auc_force_z = d$auc_force_z,
      success_z = d$success_z, age_z = d$age_z, sex_male = d$sex_male,
      rating_z = d$rating_z
    )
  )
}

summarise_effort_prior_execution <- function(fit, out) {
  effort <- draw_matrix(fit, "beta_effort_group")
  prior <- draw_matrix(fit, "beta_prior_group")
  interaction <- draw_matrix(fit, "beta_interaction_group")
  quantities <- list(
    control_required_effort = effort[, 1], MDD_required_effort = effort[, 2],
    MDD_minus_control_required_effort = effort[, 2] - effort[, 1],
    control_early_estimate = prior[, 1], MDD_early_estimate = prior[, 2],
    MDD_minus_control_early_estimate = prior[, 2] - prior[, 1],
    control_early_estimate_by_effort = interaction[, 1],
    MDD_early_estimate_by_effort = interaction[, 2],
    MDD_minus_control_early_estimate_by_effort = interaction[, 2] - interaction[, 1]
  )
  summary <- imap_dfr(quantities, ~ summarise_draws_vector(.x, .y))
  correlation <- cor(effort[, 2] - effort[, 1], interaction[, 2] - interaction[, 1])
  write_csv(summary, file.path(out, "execution_adjusted_summary.csv"))
  write_csv(tibble(
    pair = "required_effort_contrast_vs_early_estimate_by_effort_contrast",
    posterior_correlation = correlation
  ), file.path(out, "execution_adjusted_posterior_correlation.csv"))
  list(summary = summary, correlation = correlation)
}

run_effort_prior <- function(inputs) {
  out <- file.path(OUTPUT_DIR, "effort_prior")
  primary_data <- prepare_effort_prior(inputs)
  primary_fit <- sample_stan("force_prior.stan", primary_data$stan_data,
                             file.path(out, "primary"), "effort_prior_primary",
                             seed = 123L)
  primary_summary <- summarise_effort_prior(primary_fit, primary_data,
                                            file.path(out, "primary"))

  execution_data <- prepare_effort_prior_execution(primary_data)
  execution_fit <- sample_stan(
    "force_prior_execution_adjusted.stan", execution_data$stan_data,
    file.path(out, "execution_adjusted"), "effort_prior_execution_adjusted",
    seed = 20260727L
  )
  execution_summary <- summarise_effort_prior_execution(
    execution_fit, file.path(out, "execution_adjusted")
  )

  kernel20 <- NULL
  if ("early_force_mvc_20" %in% names(inputs$ep)) {
    kernel20_data <- prepare_effort_prior(inputs, "early_force_mvc_20")
    kernel20_fit <- sample_stan("force_prior.stan", kernel20_data$stan_data,
                               file.path(out, "kernel20_sensitivity"),
                               "effort_prior_kernel20", seed = 20260727L)
    kernel20_summary <- summarise_effort_prior(
      kernel20_fit, kernel20_data, file.path(out, "kernel20_sensitivity"),
      "sensitivity_20_sample"
    )
    kernel20 <- list(fit = kernel20_fit, prepared = kernel20_data,
                     summary = kernel20_summary)
  }
  list(
    primary = list(fit = primary_fit, prepared = primary_data,
                   summary = primary_summary),
    execution = list(fit = execution_fit, prepared = execution_data,
                     summary = execution_summary),
    kernel20 = kernel20
  )
}
