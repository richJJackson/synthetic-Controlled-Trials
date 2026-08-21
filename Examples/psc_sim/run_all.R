# Run all 18 scenarios with resume support (safe to stop with Ctrl+C between scenarios).

if (!exists("PSC_SIM_DIR")) {
  PSC_SIM_DIR <- normalizePath(dirname(sys.frames()[[1]]$ofile), mustWork = FALSE)
}
assign("PSC_SIM_DIR", PSC_SIM_DIR, envir = .GlobalEnv)
Sys.setenv(PSC_SIM_DIR = PSC_SIM_DIR)

source(file.path(PSC_SIM_DIR, "setup.R"))
sim_source_all(PSC_SIM_DIR)
source(file.path(PSC_SIM_DIR, "run_scenario.R"))

paths <- sim_ensure_dirs()
sim_load_packages(paths)

pop_file <- file.path(paths$data_dir, "populations.RData")
if (!file.exists(pop_file)) {
  source(file.path(PSC_SIM_DIR, "generate_populations.R"))
}

cat("Progress manifest:", sim_manifest_path(paths), "\n")
cat("Stop safely with Ctrl+C between scenarios or after a checkpoint message.\n\n")
sim_print_progress(paths)

sim_run_batch(
  populations_file = pop_file,
  resume = TRUE,
  force_restart = FALSE,
  verbose = TRUE
)
