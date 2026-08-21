# Print simulation progress without running anything.

if (!exists("PSC_SIM_DIR")) {
  PSC_SIM_DIR <- normalizePath(dirname(sys.frames()[[1]]$ofile), mustWork = FALSE)
}
assign("PSC_SIM_DIR", PSC_SIM_DIR, envir = .GlobalEnv)
Sys.setenv(PSC_SIM_DIR = PSC_SIM_DIR)

source(file.path(PSC_SIM_DIR, "setup.R"))
source(file.path(PSC_SIM_DIR, "config.R"))
source(file.path(PSC_SIM_DIR, "dgp.R"))
source(file.path(PSC_SIM_DIR, "summarise.R"))
source(file.path(PSC_SIM_DIR, "progress.R"))

sim_print_progress()
