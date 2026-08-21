# Re-apply corrected PSC MCMC to existing scenario results (PSC columns only).
#
# Usage (from psc_sim/):
#   source("run_rerun_mcmc_fix.R")
#   Rscript run_rerun_mcmc_fix.R

if (!exists("PSC_SIM_DIR")) {
  PSC_SIM_DIR <- normalizePath(getwd())
}
assign("PSC_SIM_DIR", PSC_SIM_DIR, envir = .GlobalEnv)
Sys.setenv(PSC_SIM_DIR = PSC_SIM_DIR)

source(file.path(PSC_SIM_DIR, "patch_psc_scenarios.R"))
