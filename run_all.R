source(file.path("R", "00_setup.R"))
source(file.path("R", "01_table1.R"))
source(file.path("R", "02_primary_ep.R"))
source(file.path("R", "03_ep_sensitivities.R"))
source(file.path("R", "04_ordinal_piecewise.R"))
source(file.path("R", "05_dynamic_observer.R"))
source(file.path("R", "06_effort_prior.R"))
source(file.path("R", "07_decision_making.R"))
source(file.path("R", "08_clinical_cross_task.R"))
source(file.path("R", "09_figures_tables.R"))

run_all <- function() {
  set.seed(SEED)
  inputs <- read_release_data(include_medications = TRUE)
  write_session_info()

  table1 <- run_table1(inputs)
  ep_primary <- run_primary_ep(inputs)
  ep_sensitivity <- run_ep_sensitivities(inputs)
  ordinal_piecewise <- run_ordinal_piecewise(inputs)
  observer <- run_dynamic_observer(inputs)
  effort_prior <- run_effort_prior(inputs)
  dm <- run_decision_making(inputs)
  associations <- run_clinical_cross_task(inputs, ep_primary, dm)

  build_manuscript_outputs(
    inputs = inputs,
    table1 = table1,
    ep_primary = ep_primary,
    ep_sensitivity = ep_sensitivity,
    ordinal_piecewise = ordinal_piecewise,
    observer = observer,
    effort_prior = effort_prior,
    dm = dm,
    associations = associations
  )
  invisible(TRUE)
}

if (sys.nframe() == 0L) run_all()

