# Complete scenarios missing a result file (full run, not PSC-only patch).
# Use after an interrupted full re-run or for never-started scenarios.
#
# Usage:
#   source("run_recover_missing_scenarios.R")

if (!exists("PSC_SIM_DIR")) {
  PSC_SIM_DIR <- normalizePath(getwd())
}
assign("PSC_SIM_DIR", PSC_SIM_DIR, envir = .GlobalEnv)
Sys.setenv(PSC_SIM_DIR = PSC_SIM_DIR)

source(file.path(PSC_SIM_DIR, "setup.R"))
sim_source_all(PSC_SIM_DIR)
source(file.path(PSC_SIM_DIR, "run_scenario.R"))
sim_load_packages()

paths <- sim_ensure_dirs()
pop_file <- file.path(paths$data_dir, "populations.RData")

missing <- c(
  "no_overlap_small_single_arm",
  "no_overlap_large_hybrid",
  "moderate_overlap_large_hybrid",
  "large_overlap_large_hybrid"
)

for (id in missing) {
  parsed <- parse_scenario_id(id)
  cat("\nRecovering:", id, "\n")
  sim_run_scenario(
    overlap = parsed$overlap,
    sample_size = parsed$sample_size,
    trial_type = parsed$trial_type,
    populations_file = pop_file,
    resume = TRUE,
    force_restart = FALSE,
    verbose = TRUE
  )
}

sim_combine_summaries(paths)
sim_print_progress(paths)
