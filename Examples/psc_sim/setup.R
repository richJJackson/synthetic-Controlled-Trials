# Package and PSC development setup.

sim_source_all <- function(sim_dir = NULL) {
  if (is.null(sim_dir)) {
    sim_dir <- sim_paths()$sim_dir
  }
  source(file.path(sim_dir, "config.R"), local = FALSE)
  source(file.path(sim_dir, "dgp.R"), local = FALSE)
  source(file.path(sim_dir, "psc_intervals.R"), local = FALSE)
  source(file.path(sim_dir, "psc_lik_prep.R"), local = FALSE)
  source(file.path(sim_dir, "stan_models.R"), local = FALSE)
  source(file.path(sim_dir, "summarise.R"), local = FALSE)
  source(file.path(sim_dir, "progress.R"), local = FALSE)
  source(file.path(sim_dir, "prepare_cohorts.R"), local = FALSE)
  source(file.path(sim_dir, "methods_single_arm.R"), local = FALSE)
  source(file.path(sim_dir, "methods_hybrid.R"), local = FALSE)
  source(file.path(sim_dir, "methods_psc_patch.R"), local = FALSE)
  invisible(sim_dir)
}

sim_load_packages <- function(paths = sim_paths(), load_psc_dev = TRUE) {
  suppressPackageStartupMessages({
    library(survival)
    library(flexsurv)
    library(ebal)
    library(rstan)
  })
  psc_pkg <- path.expand(paths$psc_pkg)
  if (dir.exists(psc_pkg)) {
    devtools::load_all(psc_pkg, quiet = TRUE)
  } else {
    suppressPackageStartupMessages(library(psc))
  }
  psc_rnew <- path.expand(paths$psc_rnew)
  if (load_psc_dev && file.exists(psc_rnew)) {
    source(psc_rnew, local = FALSE)
  } else if (load_psc_dev) {
    warning("PSC development file not found: ", psc_rnew)
  }
  sim_init_stan(paths, verbose = interactive())
  invisible(TRUE)
}

sim_ensure_dirs <- function(paths = sim_paths()) {
  dir.create(paths$data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
  invisible(paths)
}
