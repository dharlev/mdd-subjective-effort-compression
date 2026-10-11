ordinal_pair_table <- function(trials, successful_only = FALSE) {
  split(trials, trials$participant_id) %>%
    imap_dfr(function(d, participant_id) {
      pairs <- utils::combn(seq_len(nrow(d)), 2)
      map_dfr(seq_len(ncol(pairs)), function(k) {
        i <- pairs[1, k]; j <- pairs[2, k]
        if (d$required_effort_pct[i] == d$required_effort_pct[j]) return(NULL)
        low <- if (d$required_effort_pct[i] < d$required_effort_pct[j]) i else j
        high <- if (low == i) j else i
        if (successful_only && !(d$success[low] == 1 && d$success[high] == 1)) return(NULL)
        tibble(
          participant_id = participant_id,
          diagnosis = d$diagnosis[1],
          low_effort = d$required_effort_pct[low],
          high_effort = d$required_effort_pct[high],
          pair = paste(low_effort, high_effort, sep = "_"),
          score = case_when(
            d$rating[high] > d$rating[low] ~ 1,
            d$rating[high] < d$rating[low] ~ 0,
            TRUE ~ 0.5
          )
        )
      })
    })
}

ordinal_subject_summary <- function(pair_table, participant_covariates) {
  pair_table %>%
    group_by(participant_id, diagnosis) %>%
    summarise(ordinal_accuracy = mean(score), comparable_pairs = n(), .groups = "drop") %>%
    left_join(participant_covariates, by = c("participant_id", "diagnosis"))
}

draw_lm_coefficients <- function(model, draws = CHAINS * SAMPLING, seed = SEED) {
  set.seed(seed)
  MASS::mvrnorm(draws, coef(model), vcov(model))
}

summarise_ordinal_model <- function(subject_data, label) {
  model <- lm(ordinal_accuracy ~ diagnosis + age_z + sex_male + subject_success_rate,
              data = subject_data)
  draws <- draw_lm_coefficients(model)
  contrast <- draws[, "diagnosisMDD"]
  descriptive <- subject_data %>%
    group_by(diagnosis) %>%
    summarise(mean = mean(ordinal_accuracy), sd = sd(ordinal_accuracy), n = n(),
              .groups = "drop")
  summary <- summarise_draws_vector(contrast, "MDD_minus_control_ordinal_accuracy",
                                    scale = "proportion") %>%
    mutate(analysis = label, .before = 1)
  list(model = model, descriptive = descriptive, summary = summary)
}

add_piecewise_variables <- function(data) {
  data %>%
    mutate(
      x_linear = (required_effort_pct - 10) / 20,
      x10_30 = pmin(pmax(required_effort_pct - 10, 0), 20) / 20,
      x30_50 = pmin(pmax(required_effort_pct - 30, 0), 20) / 20,
      x50_70 = pmin(pmax(required_effort_pct - 50, 0), 20) / 20,
      x_after30 = x30_50 + x50_70,
      x_before50 = x10_30 + x30_50
    )
}

fit_piecewise_group <- function(data, diagnosis_label) {
  d <- data %>% filter(diagnosis == diagnosis_label)
  d$participant_id <- factor(d$participant_id)
  control <- nlme::lmeControl(opt = "optim", msMaxIter = 300)
  models <- list(
    one_slope = nlme::lme(rating ~ x_linear, random = ~ 1 | participant_id,
                          data = d, method = "ML", control = control),
    two_slopes_break_after_30 = nlme::lme(
      rating ~ x10_30 + x_after30, random = ~ 1 | participant_id,
      data = d, method = "ML", control = control
    ),
    two_slopes_break_after_50 = nlme::lme(
      rating ~ x_before50 + x50_70, random = ~ 1 | participant_id,
      data = d, method = "ML", control = control
    ),
    three_interval_slopes = nlme::lme(
      rating ~ x10_30 + x30_50 + x50_70, random = ~ 1 | participant_id,
      data = d, method = "ML", control = control
    )
  )
  comparison <- imap_dfr(models, function(model, name) {
    tibble(
      diagnosis = diagnosis_label, model = name, n_observations = nobs(model),
      log_likelihood = as.numeric(logLik(model)),
      parameters = attr(logLik(model), "df"), AIC = AIC(model), BIC = BIC(model)
    )
  }) %>% mutate(delta_BIC_from_best = BIC - min(BIC))
  coefficients <- imap_dfr(models, function(model, name) {
    tab <- summary(model)$tTable
    as_tibble(tab, rownames = "term") %>%
      transmute(diagnosis = diagnosis_label, model = name, term,
                estimate = Value, se = Std.Error, df = DF,
                t = `t-value`, p = `p-value`)
  })
  list(models = models, comparison = comparison, coefficients = coefficients)
}

run_ordinal_piecewise <- function(inputs) {
  out <- file.path(OUTPUT_DIR, "ordinal_piecewise")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  participants <- inputs$participants %>%
    filter(ep_primary_include == 1) %>%
    mutate(age_z = z_safe(age_years), sex_male = as.numeric(sex == "male"))
  trials <- inputs$ep %>%
    inner_join(participants %>%
                 select(participant_id, diagnosis, age_z, sex_male),
               by = "participant_id") %>%
    filter(include_primary == 1, is_practice == 0) %>%
    group_by(participant_id) %>%
    mutate(subject_success_rate = mean(success)) %>%
    ungroup()
  covariates <- trials %>%
    distinct(participant_id, diagnosis, age_z, sex_male, subject_success_rate)

  all_pairs <- ordinal_pair_table(trials, FALSE)
  successful_pairs <- ordinal_pair_table(trials, TRUE)
  all_subject <- ordinal_subject_summary(all_pairs, covariates)
  successful_subject <- ordinal_subject_summary(successful_pairs, covariates)
  ordinal_all <- summarise_ordinal_model(all_subject, "all_trials")
  ordinal_success <- summarise_ordinal_model(successful_subject, "successful_trials_only")
  ordinal_summary <- bind_rows(ordinal_all$summary, ordinal_success$summary)
  ordinal_descriptive <- bind_rows(
    ordinal_all$descriptive %>% mutate(analysis = "all_trials", .before = 1),
    ordinal_success$descriptive %>% mutate(analysis = "successful_trials_only", .before = 1)
  )
  write_csv(ordinal_summary, file.path(out, "ordinal_summary.csv"))
  write_csv(ordinal_descriptive, file.path(out, "ordinal_descriptive.csv"))

  subject_means <- trials %>%
    group_by(participant_id, diagnosis, required_effort_pct) %>%
    summarise(rating = mean(rating), .groups = "drop") %>%
    add_piecewise_variables()
  piece_control <- fit_piecewise_group(subject_means, "control")
  piece_mdd <- fit_piecewise_group(subject_means, "MDD")
  piecewise_comparison <- bind_rows(piece_control$comparison, piece_mdd$comparison)
  piecewise_coefficients <- bind_rows(piece_control$coefficients, piece_mdd$coefficients)
  write_csv(piecewise_comparison, file.path(out, "piecewise_model_comparison.csv"))
  write_csv(piecewise_coefficients, file.path(out, "piecewise_coefficients.csv"))

  list(
    ordinal_summary = ordinal_summary,
    ordinal_descriptive = ordinal_descriptive,
    ordinal_subject = list(all = all_subject, successful = successful_subject),
    piecewise_comparison = piecewise_comparison,
    piecewise_coefficients = piecewise_coefficients,
    subject_means = subject_means
  )
}

