prepare_dm <- function(inputs) {
  participants <- inputs$participants %>%
    filter(dm_model_include == 1) %>%
    arrange(source_group, participant_id) %>%
    mutate(subject_index = row_number(), group_id = as.integer(source_group))
  check_group_counts(participants, "dm_model_include", EXPECTED_COUNTS$dm,
                     "Final decision-making model")
  trials <- inputs$dm %>%
    inner_join(participants %>% select(participant_id, subject_index, group_id),
               by = "participant_id") %>%
    filter(is_practice == 0, choice_option2 %in% c(0, 1)) %>%
    arrange(subject_index, trial)
  reward_max <- max(c(trials$reward1, trials$reward2), na.rm = TRUE)
  trials <- trials %>%
    mutate(
      reward1_norm = reward1 / reward_max,
      reward2_norm = reward2 / reward_max,
      effort1_norm = effort1_pct / 100,
      effort2_norm = effort2_pct / 100
    )
  stan_data <- list(
    N_subj = nrow(participants), N_obs = nrow(trials),
    G = length(GROUP_LEVELS), subj = trials$subject_index,
    group_id = participants$group_id, y = as.integer(trials$choice_option2),
    rew1 = trials$reward1_norm, eff1 = trials$effort1_norm,
    rew2 = trials$reward2_norm, eff2 = trials$effort2_norm
  )
  list(participants = participants, trials = trials, reward_max = reward_max,
       stan_data = stan_data)
}

summarise_dm <- function(fit, prepared, out, prefix = "dm") {
  p <- prepared$participants
  weights <- control_weights(p, "dm_model_include")
  mu_reward <- draw_matrix(fit, "mu_rew")
  mu_effort <- draw_matrix(fit, "mu_eff")
  reward <- draw_matrix(fit, "rewSens")
  effort <- draw_matrix(fit, "effSens")
  control_reward <- weighted_control_draw(mu_reward, weights)
  control_effort <- weighted_control_draw(mu_effort, weights)
  quantities <- list(
    control_reward_sensitivity = control_reward,
    MDD_reward_sensitivity = mu_reward[, 3],
    MDD_minus_control_reward_sensitivity = mu_reward[, 3] - control_reward,
    control_effort_sensitivity = control_effort,
    MDD_effort_sensitivity = mu_effort[, 3],
    MDD_minus_control_effort_sensitivity = mu_effort[, 3] - control_effort
  )
  summary <- imap_dfr(quantities, ~ summarise_draws_vector(.x, .y, "model scale"))
  participant_means <- p %>%
    transmute(
      participant_id, diagnosis, source_group,
      reward_sensitivity_posterior_mean = colMeans(reward),
      effort_sensitivity_posterior_mean = colMeans(effort)
    )
  ppc_draws <- as.numeric(draw_matrix(fit, "ybar_rep"))
  ppc <- tibble(
    observed_option2_rate = mean(prepared$stan_data$y),
    posterior_predicted_option2_rate = mean(ppc_draws),
    q025 = quantile(ppc_draws, 0.025),
    q975 = quantile(ppc_draws, 0.975)
  )
  write_csv(summary, file.path(out, paste0(prefix, "_group_summary.csv")))
  write_csv(participant_means, file.path(out, paste0(prefix, "_participant_means.csv")))
  write_csv(ppc, file.path(out, paste0(prefix, "_posterior_predictive.csv")))
  list(summary = summary, participant_means = participant_means,
       reward_draws = reward, effort_draws = effort, ppc = ppc,
       ppc_draws = ppc_draws)
}

run_dm_recovery <- function(empirical_fit, prepared, out) {
  array_draws <- empirical_fit$draws(c("rewSens", "effSens"), format = "draws_array")
  dims <- dim(array_draws)
  flattened <- do.call(rbind, lapply(seq_len(dims[2]), function(chain) {
    array_draws[, chain, , drop = FALSE][, 1, ]
  }))
  reward <- flattened[, grepl("^rewSens\\[", colnames(flattened)), drop = FALSE]
  effort <- flattened[, grepl("^effSens\\[", colnames(flattened)), drop = FALSE]
  set.seed(123L)
  generating_draw <- sample.int(nrow(flattened), 1)
  true_reward <- reward[generating_draw, ]
  true_effort <- effort[generating_draw, ]
  d <- prepared$stan_data
  probability <- vapply(seq_len(d$N_obs), function(i) {
    s <- d$subj[i]
    value1 <- true_reward[s] * d$rew1[i] - true_effort[s] * d$eff1[i]^2
    value2 <- true_reward[s] * d$rew2[i] - true_effort[s] * d$eff2[i]^2
    plogis(value2 - value1)
  }, numeric(1))
  set.seed(123L)
  invisible(sample.int(nrow(flattened), 1))
  recovery_data <- d
  recovery_data$y <- rbinom(d$N_obs, 1, probability)
  recovery_fit <- sample_stan("decision_making.stan", recovery_data, out,
                              "dm_parameter_recovery", seed = 456L)
  recovered_reward <- colMeans(draw_matrix(recovery_fit, "rewSens"))
  recovered_effort <- colMeans(draw_matrix(recovery_fit, "effSens"))
  recovery <- prepared$participants %>%
    transmute(
      participant_id, diagnosis, source_group,
      true_reward_sensitivity = true_reward,
      recovered_reward_sensitivity = recovered_reward,
      true_effort_sensitivity = true_effort,
      recovered_effort_sensitivity = recovered_effort
    )
  recovery_summary <- tibble(
    parameter = c("reward_sensitivity", "effort_sensitivity"),
    correlation = c(cor(true_reward, recovered_reward),
                    cor(true_effort, recovered_effort)),
    median_absolute_log_error = c(
      median(abs(log(recovered_reward) - log(true_reward))),
      median(abs(log(recovered_effort) - log(true_effort)))
    )
  )
  write_csv(recovery, file.path(out, "dm_recovery_participants.csv"))
  write_csv(recovery_summary, file.path(out, "dm_recovery_summary.csv"))
  list(fit = recovery_fit, participant = recovery, summary = recovery_summary)
}

run_decision_making <- function(inputs) {
  out <- file.path(OUTPUT_DIR, "decision_making")
  prepared <- prepare_dm(inputs)
  fit <- sample_stan("decision_making.stan", prepared$stan_data, out,
                     "decision_making", seed = 123L)
  summary <- summarise_dm(fit, prepared, out)
  recovery <- NULL
  if (tolower(Sys.getenv("RUN_DM_RECOVERY", "true")) == "true") {
    recovery <- run_dm_recovery(fit, prepared, file.path(out, "recovery"))
  }
  list(fit = fit, prepared = prepared, summary = summary, recovery = recovery)
}
