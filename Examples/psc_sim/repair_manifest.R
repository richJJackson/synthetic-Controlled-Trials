# Rebuild a corrupted progress manifest from files on disk.
# Use if show_progress.R errors or scenario names look wrong.

if (!exists("PSC_SIM_DIR")) {
  PSC_SIM_DIR <- normalizePath(dirname(sys.frames()[[1]]$ofile), mustWork = FALSE)
}
Sys.setenv(PSC_SIM_DIR = PSC_SIM_DIR)

source(file.path(PSC_SIM_DIR, "setup.R"))
source(file.path(PSC_SIM_DIR, "config.R"))
source(file.path(PSC_SIM_DIR, "dgp.R"))
source(file.path(PSC_SIM_DIR, "summarise.R"))
source(file.path(PSC_SIM_DIR, "progress.R"))

paths <- sim_ensure_dirs()
manifest_file <- sim_manifest_path(paths)

if (file.exists(manifest_file)) {
  backup <- paste0(manifest_file, ".bak_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  file.copy(manifest_file, backup)
  cat("Backed up old manifest to:\n ", backup, "\n", sep = "")
}

manifest <- sim_init_manifest(paths = paths)
manifest <- sim_sync_manifest_from_disk(manifest, paths)
cat("Manifest rebuilt:", manifest_file, "\n\n")
sim_print_progress(paths, sync = FALSE)
