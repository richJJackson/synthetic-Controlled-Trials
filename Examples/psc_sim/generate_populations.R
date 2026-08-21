# Generate super-populations and save to data/populations.RData.

if (!exists("PSC_SIM_DIR")) {
  PSC_SIM_DIR <- normalizePath(dirname(sys.frames()[[1]]$ofile), mustWork = FALSE)
}
assign("PSC_SIM_DIR", PSC_SIM_DIR, envir = .GlobalEnv)
Sys.setenv(PSC_SIM_DIR = PSC_SIM_DIR)

source(file.path(PSC_SIM_DIR, "setup.R"))
sim_source_all(PSC_SIM_DIR)
paths <- sim_ensure_dirs()

populations <- sim_generate_all_populations()
out_file <- file.path(paths$data_dir, "populations.RData")
save(populations, file = out_file)

cat("Saved super-populations:", out_file, "\n")
cat("  control: n =", nrow(populations$control), "\n")
for (ov in SIM_OVERLAP_LEVELS) {
  nm <- trial_population_name(ov)
  cat(" ", nm, ": n =", nrow(populations[[nm]]), "\n")
}
