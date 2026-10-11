cran_packages <- c(
  "posterior", "dplyr", "tidyr", "readr", "tibble", "purrr",
  "ggplot2", "patchwork", "nlme", "MASS", "loo", "scales"
)
install.packages(setdiff(cran_packages, rownames(installed.packages())),
                 repos = "https://cloud.r-project.org")

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages(
    "cmdstanr",
    repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org")
  )
}

cmdstan_available <- tryCatch({
  cmdstanr::cmdstan_version()
  TRUE
}, error = function(e) FALSE)
if (!cmdstan_available) {
  cmdstanr::install_cmdstan()
}
