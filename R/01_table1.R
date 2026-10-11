cohen_d_equal_variance <- function(x_mdd, x_control) {
  n1 <- sum(is.finite(x_mdd)); n0 <- sum(is.finite(x_control))
  pooled <- sqrt(((n1 - 1) * var(x_mdd, na.rm = TRUE) +
                    (n0 - 1) * var(x_control, na.rm = TRUE)) / (n1 + n0 - 2))
  (mean(x_mdd, na.rm = TRUE) - mean(x_control, na.rm = TRUE)) / pooled
}

numeric_table1_row <- function(data, variable, label) {
  control <- data[[variable]][data$diagnosis == "control"]
  mdd <- data[[variable]][data$diagnosis == "MDD"]
  test <- t.test(mdd, control, var.equal = TRUE)
  tibble(
    characteristic = label,
    control = sprintf("%.2f (%.2f)", mean(control, na.rm = TRUE), sd(control, na.rm = TRUE)),
    MDD = sprintf("%.2f (%.2f)", mean(mdd, na.rm = TRUE), sd(mdd, na.rm = TRUE)),
    statistic = unname(test$statistic),
    df = unname(test$parameter),
    p = test$p.value,
    effect = cohen_d_equal_variance(mdd, control),
    effect_name = "Cohen d (MDD - control)"
  )
}

sex_table1_row <- function(data) {
  tab <- table(data$diagnosis, data$sex)
  test <- chisq.test(tab, correct = FALSE)
  n_control <- sum(tab["control", ])
  n_mdd <- sum(tab["MDD", ])
  phi <- sqrt(unname(test$statistic) / sum(tab))
  tibble(
    characteristic = "Female sex, n (%)",
    control = sprintf("%d (%.1f%%)", tab["control", "female"],
                      100 * tab["control", "female"] / n_control),
    MDD = sprintf("%d (%.1f%%)", tab["MDD", "female"],
                  100 * tab["MDD", "female"] / n_mdd),
    statistic = unname(test$statistic),
    df = unname(test$parameter),
    p = test$p.value,
    effect = phi,
    effect_name = "Phi"
  )
}

run_table1 <- function(inputs) {
  out <- file.path(OUTPUT_DIR, "tables")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  p <- inputs$participants %>% filter(table1_include == 1)
  check_group_counts(p, "table1_include", EXPECTED_COUNTS$table1, "Table 1")

  rows <- list(
    numeric_table1_row(p, "age_years", "Age, years"),
    sex_table1_row(p),
    numeric_table1_row(p, "phq9_total", "PHQ-9"),
    numeric_table1_row(p, "aes_total", "AES")
  )
  table1 <- bind_rows(rows)
  write_csv(table1, file.path(out, "table1.csv"))

  flow <- tibble(
    sample = c("Contributed to any analysis", "Primary EP", "Final DM model"),
    controls = c(
      sum(p$diagnosis == "control"),
      sum(inputs$participants$ep_primary_include == 1 & inputs$participants$diagnosis == "control"),
      sum(inputs$participants$dm_model_include == 1 & inputs$participants$diagnosis == "control")
    ),
    MDD = c(
      sum(p$diagnosis == "MDD"),
      sum(inputs$participants$ep_primary_include == 1 & inputs$participants$diagnosis == "MDD"),
      sum(inputs$participants$dm_model_include == 1 & inputs$participants$diagnosis == "MDD")
    )
  )
  write_csv(flow, file.path(out, "participant_flow.csv"))
  list(table = table1, flow = flow)
}

