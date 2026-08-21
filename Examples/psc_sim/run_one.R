# Run a single scenario from the command line or interactively.
#
# Examples:
#   source("run_one.R"); run_one_scenario("no_overlap", "small", "single_arm")
#   run_one_scenario("no_overlap", "small", "single_arm", force_restart = TRUE)
#   Rscript run_one.R no_overlap small hybrid

if (!exists("PSC_SIM_DIR")) {
  args <- commandArgs(trailingOnly = TRUE)
  PSC_SIM_DIR <- normalizePath(
    if (length(args) >= 4) args[4] else dirname(sys.frames()[[1]]$ofile),
    mustWork = FALSE
  )
}
assign("PSC_SIM_DIR", PSC_SIM_DIR, envir = .GlobalEnv)
Sys.setenv(PSC_SIM_DIR = PSC_SIM_DIR)

source(file.path(PSC_SIM_DIR, "setup.R"))
sim_source_all(PSC_SIM_DIR)
source(file.path(PSC_SIM_DIR, "run_scenario.R"))
sim_load_packages()

run_one_scenario <- function(
    overlap = "no_overlap",
    sample_size = "small",
    trial_type = "single_arm",
    nsim = NULL,
  ...) {
  sim_run_scenario(
    overlap = overlap,
    sample_size = sample_size,
    trial_type = trial_type,
    nsim = nsim,
    verbose = TRUE,
    ...
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  run_one_scenario(args[1], args[2], args[3], if (length(args) >= 4) as.integer(args[4]) else NULL)
}
