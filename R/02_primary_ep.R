prepare_primary_ep <- function(inputs, ceiling_threshold = Inf) {
  participants <- inputs$participants %>%
    filter(ep_primary_include == 1) %>%
    arrange(source_group, participant_id) %>%
    mutate(subject_index = row_number(), group_id = as.integer(source_group))
  check_group_counts(participants, "ep_primary_include", EXPECTED_COUNTS$ep,
                     "Primary EP")

  trials <- inputs$ep %>%
    inner_join(participants %>% select(participant_id, subject_index, group_id),
               by = "participant_id") %>%
    filter(include_primary == 1, is_practice == 0,
           is.finite(required_effort_pct), is.finite(rating), is.finite(trial),
           rating <= ceiling_threshold) %>%
    arrange(subject_index, trial)

  scaling <- tibble(
    variable = c("rating", "required_effort_pct", "trial"),
    mean = c(mean(trials$rating), mean(trials$required_effort_pct), mean(trials$trial)),
    sd = c(sd(trials$rating), sd(trials$required_effort_pct), sd(trials$trial))
  )
  trials <- trials %>%
    mutate(
      rating_z = z_safe(rating),
      effort_z = z_safe(required_effort_pct),
      trial_z = z_safe(trial)
    )
  stan_data <- list(
    N_subj = nrow(participants),
    N_obs = nrow(trials),
    N_group = length(GROUP_LEVELS),
    subj = trials$subject_index,
    subj_group = participants$group_id,
    rating_z = trials$rating_z,
    effort_z = trials$effort_z,
    trial_z = trials$trial_z
  )
  list(participants = participants, trials = trials, scaling = scaling,
       stan_data = stan_data)
}

summarise_primary_ep <- function(fit, prepared, output_dir, prefix = "primary") {
  participants <- prepared$participants
  alpha <- draw_matrix(fit, "alpha")
  beta <- draw_matrix(fit, "beta")
  mu_alpha <- draw_matrix(fit, "mu_alpha_group")
  mu_beta <- draw_matrix(fit, "mu_beta_group")

  idx_control <- which(participants$diagnosis == "control")
  idx_mdd <- which(participants$diagnosis == "MDD")
  weights <- control_weights(participants, "ep_primary_include")

  finite_draws <- list(
    control_intercept = rowMeans(alpha[, idx_control, drop = FALSE]),
    MDD_intercept = rowMeans(alpha[, idx_mdd, drop = FALSE]),
    control_discriminability = rowMeans(beta[, idx_control, drop = FALSE]),
    MDD_discriminability = rowMeans(beta[, idx_mdd, drop = FALSE])
  )
  finite_draws$MDD_minus_control_intercept <-
    finite_draws$MDD_intercept - finite_draws$control_intercept
  finite_draws$MDD_minus_control_discriminability <-
    finite_draws$MDD_discriminability - finite_draws$control_discriminability

  population_draws <- list(
    control_intercept = weighted_control_draw(mu_alpha, weights),
    MDD_intercept = mu_alpha[, 3],
    control_discriminability = weighted_control_draw(mu_beta, weights),
    MDD_discriminability = mu_beta[, 3]
  )
  population_draws$MDD_minus_control_intercept <-
    population_draws$MDD_intercept - population_draws$control_intercept
  population_draws$MDD_minus_control_discriminability <-
    population_draws$MDD_discriminability - population_draws$control_discriminability

  finite_summary <- imap_dfr(
    finite_draws,
    ~ summarise_draws_vector(.x, .y) %>%
      mutate(estimand = "finite_sample_participant_average", .before = 1)
  )
  population_summary <- imap_dfr(
    population_draws,
    ~ summarise_draws_vector(.x, .y) %>%
      mutate(estimand = "hierarchical_population_mean", .before = 1)
  )
  summary <- bind_rows(finite_summary, population_summary)

  scaling <- setNames(prepared$scaling$sd, prepared$scaling$variable)
  raw_scale <- summary %>%
    mutate(
      mean_raw = case_when(
        grepl("discriminability", term) ~ mean * scaling[["rating"]] /
          scaling[["required_effort_pct"]] * 10,
        TRUE ~ mean * scaling[["rating"]]
      ),
      q025_raw = case_when(
        grepl("discriminability", term) ~ q025 * scaling[["rating"]] /
          scaling[["required_effort_pct"]] * 10,
        TRUE ~ q025 * scaling[["rating"]]
      ),
      q975_raw = case_when(
        grepl("discriminability", term) ~ q975 * scaling[["rating"]] /
          scaling[["required_effort_pct"]] * 10,
        TRUE ~ q975 * scaling[["rating"]]
      ),
      raw_unit = if_else(grepl("discriminability", term),
                         "rating points per 10% MVC", "rating points")
    )

  participant_means <- prepared$participants %>%
    transmute(
      participant_id, diagnosis, source_group,
      intercept_posterior_mean = colMeans(alpha),
      discriminability_posterior_mean = colMeans(beta)
    )
  write_csv(summary, file.path(output_dir, paste0(prefix, "_estimands.csv")))
  write_csv(raw_scale, file.path(output_dir, paste0(prefix, "_estimands_with_raw_scale.csv")))
  write_csv(participant_means,
            file.path(output_dir, paste0(prefix, "_participant_posterior_means.csv")))
  write_csv(prepared$scaling, file.path(output_dir, paste0(prefix, "_scaling.csv")))
  list(summary = summary, raw_scale = raw_scale, participant_means = participant_means,
       alpha_draws = alpha, beta_draws = beta,
       finite_draws = finite_draws, population_draws = population_draws)
}

run_primary_ep <- function(inputs) {
  out <- file.path(OUTPUT_DIR, "ep_primary")
  prepared <- prepare_primary_ep(inputs)
  fit <- sample_stan("ep_primary.stan", prepared$stan_data, out, "ep_primary",
                     seed = 20260601L)
  summary <- summarise_primary_ep(fit, prepared, out)
  list(fit = fit, prepared = prepared, summary = summary)
}
