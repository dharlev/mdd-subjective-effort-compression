options(stringsAsFactors = FALSE, scipen = 999)

required_packages <- c(
  "cmdstanr", "posterior", "dplyr", "tidyr", "readr", "tibble",
  "purrr", "ggplot2", "patchwork", "nlme", "MASS", "loo", "scales"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                              logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Install the required R packages before running the pipeline: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(posterior)
  library(ggplot2)
})

REPO_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(REPO_ROOT, "run_all.R"))) {
  stop("Run the pipeline from the repository root.", call. = FALSE)
}

DATA_DIR <- file.path(REPO_ROOT, "data")
OUTPUT_DIR <- file.path(REPO_ROOT, "outputs")
MODEL_DIR <- file.path(REPO_ROOT, "models")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

SEED <- 20260723L
CHAINS <- 4L
WARMUP <- 1000L
SAMPLING <- 4000L
PARALLEL_CHAINS <- min(CHAINS, parallel::detectCores(logical = FALSE))

GROUP_LEVELS <- c("control_cohort_1", "control_cohort_2", "MDD")
DIAGNOSIS_LEVELS <- c("control", "MDD")
EXPECTED_COUNTS <- list(
  table1 = c(control = 94L, MDD = 50L),
  ep = c(control = 93L, MDD = 49L),
  dm = c(control = 93L, MDD = 41L)
)

`%||%` <- function(x, y) if (is.null(x)) y else x

z_safe <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) stop("Cannot standardise a constant variable.")
  (x - mean(x, na.rm = TRUE)) / s
}

assert_columns <- function(data, columns, label) {
  absent <- setdiff(columns, names(data))
  if (length(absent)) {
    stop(label, " is missing: ", paste(absent, collapse = ", "), call. = FALSE)
  }
}

assert_anonymous_ids <- function(x) {
  bad <- unique(x[is.na(x) | !grepl("^P[0-9]{4}$", x)])
  if (length(bad)) {
    stop("participant_id must use anonymous P#### values only.", call. = FALSE)
  }
}

read_release_data <- function(include_medications = FALSE) {
  required <- c("participants.csv", "ep_trials.csv", "dm_trials.csv")
  absent <- required[!file.exists(file.path(DATA_DIR, required))]
  if (length(absent)) {
    stop("Missing input file(s): ", paste(absent, collapse = ", "),
         ". See DATA_SCHEMA.md.", call. = FALSE)
  }

  participants <- read_csv(file.path(DATA_DIR, "participants.csv"), show_col_types = FALSE)
  ep <- read_csv(file.path(DATA_DIR, "ep_trials.csv"), show_col_types = FALSE)
  dm <- read_csv(file.path(DATA_DIR, "dm_trials.csv"), show_col_types = FALSE)

  assert_columns(participants, c(
    "participant_id", "diagnosis", "control_cohort", "age_years", "sex",
    "phq9_total", "aes_total", "ami_total", "table1_include",
    "ep_primary_include", "dm_task_complete", "dm_model_include"
  ), "participants.csv")
  assert_columns(ep, c(
    "participant_id", "trial", "is_practice", "required_effort_pct",
    "rating", "success", "force_auc", "peak_force", "mean_force",
    "early_force_mvc", "early_force_time_samples", "include_primary",
    "include_early_force"
  ), "ep_trials.csv")
  assert_columns(dm, c(
    "participant_id", "trial", "is_practice", "is_catch", "choice_option2",
    "reward1", "effort1_pct", "reward2", "effort2_pct", "success"
  ), "dm_trials.csv")

  assert_anonymous_ids(participants$participant_id)
  assert_anonymous_ids(ep$participant_id)
  assert_anonymous_ids(dm$participant_id)
  if (anyDuplicated(participants$participant_id)) {
    stop("participants.csv must contain one row per participant.", call. = FALSE)
  }
  if (!all(ep$participant_id %in% participants$participant_id) ||
      !all(dm$participant_id %in% participants$participant_id)) {
    stop("Trial tables contain participant IDs absent from participants.csv.", call. = FALSE)
  }

  participants <- participants %>%
    mutate(
      diagnosis = factor(diagnosis, DIAGNOSIS_LEVELS),
      sex = factor(sex, c("female", "male")),
      source_group = case_when(
        diagnosis == "MDD" ~ "MDD",
        control_cohort == "cohort_1" ~ "control_cohort_1",
        control_cohort == "cohort_2" ~ "control_cohort_2",
        TRUE ~ NA_character_
      ),
      source_group = factor(source_group, GROUP_LEVELS)
    )
  if (anyNA(participants$diagnosis) || anyNA(participants$sex) ||
      anyNA(participants$source_group)) {
    stop("Invalid diagnosis, sex, or control_cohort coding.", call. = FALSE)
  }

  medication_path <- file.path(DATA_DIR, "medications.csv")
  medications <- NULL
  if (include_medications && file.exists(medication_path)) {
    medications <- read_csv(medication_path, show_col_types = FALSE)
    assert_columns(medications, c(
      "participant_id", "medication_class", "medication_name", "dose_value",
      "dose_unit", "prn"
    ), "medications.csv")
    assert_anonymous_ids(medications$participant_id)
  }
  list(participants = participants, ep = ep, dm = dm, medications = medications)
}

check_group_counts <- function(participants, flag, expected, label) {
  observed <- participants %>%
    filter(.data[[flag]] == 1) %>%
    count(diagnosis) %>%
    tidyr::complete(diagnosis = factor(DIAGNOSIS_LEVELS, DIAGNOSIS_LEVELS),
                    fill = list(n = 0L))
  got <- setNames(observed$n, as.character(observed$diagnosis))
  if (!identical(as.integer(got[names(expected)]), as.integer(expected))) {
    stop(label, " counts differ from the final manuscript sample. Observed: ",
         paste(names(got), got, sep = "=", collapse = ", "), call. = FALSE)
  }
  invisible(got)
}

posterior_pd <- function(x) max(mean(x > 0), mean(x < 0))

summarise_draws_vector <- function(x, term, scale = "standardised") {
  tibble(
    term = term,
    scale = scale,
    mean = mean(x),
    sd = sd(x),
    q025 = unname(quantile(x, 0.025)),
    q975 = unname(quantile(x, 0.975)),
    pd = posterior_pd(x)
  )
}

draw_matrix <- function(fit, variable) {
  as_draws_matrix(fit$draws(variable))
}

fit_diagnostics <- function(fit, label) {
  s <- fit$summary()
  sampler <- fit$sampler_diagnostics(format = "draws_matrix")
  divergences <- sum(sampler[, , "divergent__"])
  treedepth_limit <- fit$metadata()$max_treedepth %||% 10
  treedepth_warnings <- sum(sampler[, , "treedepth__"] >= treedepth_limit)
  tibble(
    fit = label,
    chains = fit$metadata()$chains,
    warmup_per_chain = fit$metadata()$iter_warmup,
    sampling_per_chain = fit$metadata()$iter_sampling,
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_bulk_ess = min(s$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(s$ess_tail, na.rm = TRUE),
    divergences = divergences,
    treedepth_warnings = treedepth_warnings
  )
}

sample_stan <- function(model_file, data, output_dir, fit_name,
                        seed = SEED, sampling = SAMPLING) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  model <- cmdstanr::cmdstan_model(file.path(MODEL_DIR, model_file))
  fit <- model$sample(
    data = data,
    seed = seed,
    chains = CHAINS,
    parallel_chains = PARALLEL_CHAINS,
    iter_warmup = WARMUP,
    iter_sampling = sampling,
    adapt_delta = 0.95,
    max_treedepth = 12,
    refresh = 200,
    output_dir = output_dir,
    save_warmup = FALSE
  )
  saveRDS(fit, file.path(output_dir, paste0(fit_name, ".rds")))
  write_csv(fit_diagnostics(fit, fit_name),
            file.path(output_dir, paste0(fit_name, "_diagnostics.csv")))
  fit
}

control_weights <- function(participants, inclusion_flag) {
  n <- participants %>%
    filter(.data[[inclusion_flag]] == 1, diagnosis == "control") %>%
    count(source_group) %>%
    filter(source_group %in% c("control_cohort_1", "control_cohort_2"))
  w <- setNames(n$n / sum(n$n), as.character(n$source_group))
  c(control_cohort_1 = w[["control_cohort_1"]],
    control_cohort_2 = w[["control_cohort_2"]])
}

weighted_control_draw <- function(group_draws, weights) {
  group_draws[, 1] * weights[["control_cohort_1"]] +
    group_draws[, 2] * weights[["control_cohort_2"]]
}

write_session_info <- function() {
  versions <- tibble(
    package = c("R", required_packages, "CmdStan"),
    version = c(
      as.character(getRversion()),
      vapply(required_packages, function(x) as.character(packageVersion(x)), character(1)),
      tryCatch(cmdstanr::cmdstan_version(), error = function(e) NA_character_)
    )
  )
  write_csv(versions, file.path(OUTPUT_DIR, "software_versions.csv"))
  capture.output(sessionInfo(), file = file.path(OUTPUT_DIR, "sessionInfo.txt"))
}
