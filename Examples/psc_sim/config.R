# Simulation study configuration
# Source after setwd() to psc_sim/ or use sim_paths().

sim_paths <- function() {
  sim_dir <- Sys.getenv("PSC_SIM_DIR", unset = "")
  if (!nzchar(sim_dir)) {
    sim_dir <- "~/Documents/GitHub/synthetic-Controlled-Trials/Examples/psc_sim"
  }
  sim_dir <- path.expand(sim_dir)
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)
  }
  sim_dir <- normalizePath(sim_dir, mustWork = TRUE)
  list(
    sim_dir = sim_dir,
    data_dir = file.path(sim_dir, "data"),
    results_dir = file.path(sim_dir, "results"),
    psc_pkg = Sys.getenv("PSC_PKG", "~/Documents/GitHub/psc"),
    psc_rnew = Sys.getenv(
      "PSC_RNEW",
      "~/Documents/GitHub/pscDevelop/ProximityWeights/methods/load_psc_dev.R"
    )
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

SIM_OVERLAP_LEVELS <- c("no_overlap", "moderate_overlap", "large_overlap")

SIM_SAMPLE_SIZES <- list(
  small = list(n_cont = 300L, n_trt = 50L),
  medium = list(n_cont = 400L, n_trt = 75L),
  large = list(n_cont = 500L, n_trt = 100L)
)

SIM_TRIAL_TYPES <- c("single_arm", "hybrid")

SIM_DGP <- list(
  n_pop = 5000L,
  seed_pop = 2026L,
  gam_0 = 0.05,
  gam_1 = log(0.6),
  gam_2 = log(0.8),
  gam_3 = log(0.6),
  gam_4 = log(1.2),
  gam_5 = log(0.9),
  beta = log(0.7),
  censor_rate = 0.02,
  admin_censor_years = 6
)

# Covariate generators for trial super-populations.
# Control super-population always uses large_overlap covariate margins.
# Trial populations vary by overlap scenario:
#   no_overlap       — original pscDesign_ex.R treated DGP (poor transportability)
#   moderate_overlap — intermediate shift (new)
#   large_overlap    — same covariate distribution as control (good transportability)
SIM_OVERLAP_SPECS <- list(
  no_overlap = list(
    X_1 = c(mean = 1.5, sd = 1),
    X_2 = c(mean = 0.5, sd = 1),
    X_3_p = 0.6,
    X_4_p = 0.7,
    X_5_p = 0.8
  ),
  moderate_overlap = list(
    X_1 = c(mean = 0.75, sd = 1),
    X_2 = c(mean = 0.25, sd = 1),
    X_3_p = 0.45,
    X_4_p = 0.55,
    X_5_p = 0.65
  ),
  large_overlap = list(
    X_1 = c(mean = 0, sd = 1),
    X_2 = c(mean = 0, sd = 1),
    X_3_p = 0.3,
    X_4_p = 0.4,
    X_5_p = 0.5
  )
)

SIM_DEFAULTS <- list(
  nsim = 1000L,
  nsim_hybrid = 1000L,
  stan_iter = 3000L,
  stan_chains = 2L,
  stan_seed = 21319L,
  psc_nchain = 1L,
  psc_nsim = 5000L,
  verbose = FALSE,
  save_results = TRUE,
  checkpoint_every = 25L,
  parallel_cores = 4L,
  max_parallel_cores = 4L
)

# Columns stored per method per replicate: estimate, SE, interval lower, upper.
# PSC methods use 95% HPD bounds in lo/hi; other methods use Wald limits.
SIM_RES_COLS <- c("est", "se", "lo", "hi")
SIM_RES_NCOL <- length(SIM_RES_COLS)

SA_METHOD_NAMES <- c(
  "Pooled (unadj)", "Pooled (full adj)", "Pooled (partial adj)",
  "PSC full", "PSC partial",
  "SC full", "SC partial",
  "Bayes hist (info)", "Bayes hist (case-weighted)",
  "PSC full (pw)", "PSC partial (pw)"
)

HYBRID_METHOD_NAMES <- c(
  "RCT (unadj)", "RCT (adj)",
  "Pooled (unadj)", "Pooled (adj)",
  "PSC full (combined)", "PSC partial (combined)",
  "SC full", "SC partial",
  "PSC full (pw, combined)", "PSC partial (pw, combined)",
  "Bayes vague", "Bayes informative", "Bayes commensurate",
  "Bayes case-weighted", "Bayes partial case-weighted"
)

scenario_grid <- function() {
  expand.grid(
    overlap = SIM_OVERLAP_LEVELS,
    sample_size = names(SIM_SAMPLE_SIZES),
    trial_type = SIM_TRIAL_TYPES,
    stringsAsFactors = FALSE
  )
}

scenario_id <- function(overlap, sample_size, trial_type) {
  paste(overlap, sample_size, trial_type, sep = "_")
}

parse_scenario_id <- function(id) {
  for (tt in SIM_TRIAL_TYPES) {
    suffix <- paste0("_", tt)
    if (endsWith(id, suffix)) {
      rest <- substr(id, 1L, nchar(id) - nchar(suffix))
      for (ss in names(SIM_SAMPLE_SIZES)) {
        suffix2 <- paste0("_", ss)
        if (endsWith(rest, suffix2)) {
          overlap <- substr(rest, 1L, nchar(rest) - nchar(suffix2))
          return(list(
            overlap = overlap,
            sample_size = ss,
            trial_type = tt
          ))
        }
      }
    }
  }
  stop("Unrecognised scenario id: ", id)
}

print_scenario_grid <- function() {
  grid <- scenario_grid()
  grid$id <- mapply(scenario_id, grid$overlap, grid$sample_size, grid$trial_type)
  grid$n_cont <- vapply(grid$sample_size, function(s) SIM_SAMPLE_SIZES[[s]]$n_cont, integer(1))
  grid$n_trt <- vapply(grid$sample_size, function(s) SIM_SAMPLE_SIZES[[s]]$n_trt, integer(1))
  print(grid[, c("id", "overlap", "sample_size", "trial_type", "n_cont", "n_trt")], row.names = FALSE)
  invisible(grid)
}
