# Validate PSC MCMC against Stan (fixed CFM) on one simulation replicate.
#
# Usage (from psc_sim/):
#   source("validate_psc_mcmc.R")
#   validate_psc_mcmc()

validate_psc_mcmc <- function(
    overlap = "no_overlap",
    sample_size = "small",
    seed = 4242L,
    nsim = 5000L,
    burn = 500L,
    thin = 2L,
    stan_iter = 4000L,
    stan_chains = 4L,
    paths = NULL,
    save_results = TRUE) {

  sim_dir <- Sys.getenv("PSC_SIM_DIR", unset = normalizePath("."))
  if (!exists("sim_paths", mode = "function")) {
    source(file.path(sim_dir, "setup.R"), local = FALSE)
    sim_source_all(sim_dir)
  }
  if (is.null(paths)) {
    paths <- sim_paths()
  }

  source(file.path(paths$sim_dir, "psc_lik_prep.R"), local = FALSE)
  sim_load_packages()

  pop_file <- file.path(paths$data_dir, "populations.RData")
  if (!file.exists(pop_file)) {
    source(file.path(paths$sim_dir, "generate_populations.R"))
  }
  load(pop_file)
  sizes <- SIM_SAMPLE_SIZES[[sample_size]]
  cohorts <- sim_sample_cohorts(
    populations, overlap, sizes$n_cont, sizes$n_trt, seed = seed
  )

  prep <- sim_prepare_single_arm(cohorts$contDat, cohorts$trtDat)
  aligned <- sim_factor_align(prep$contDat, prep$trtDat)
  models <- sim_fit_cfm(aligned$contDat)
  trtDat <- aligned$trtDat

  pscOb <- pscData(models$cfmw, trtDat)
  pscOb <- init(pscOb)
  lik_prep <- psc_lik_prep_fixed(pscOb, co = pscOb$co)

  stan_data <- list(
    N = nrow(lik_prep),
    H0 = lik_prep$H0,
    h0 = lik_prep$h0,
    lp = lik_prep$lp,
    status = lik_prep$status
  )

  stan_fit <- sim_stan_fit(
    "psc_fixed",
    stan_data,
    defaults = modifyList(SIM_DEFAULTS, list(
      stan_iter = stan_iter,
      stan_chains = stan_chains,
      stan_seed = seed
    ))
  )
  stan_drs <- posterior::as_draws_df(stan_fit)
  stan_beta <- stan_drs$beta

  psc_fixed <- pscfit(
    models$cfmw,
    trtDat,
    nsim = nsim,
    nchain = 1L,
    burn = burn,
    thin = thin,
    fixed_cfm = TRUE,
    tune_burn = min(500L, nsim %/% 2L)
  )
  psc_full <- pscfit(
    models$cfmw,
    trtDat,
    nsim = nsim,
    nchain = 1L,
    burn = burn,
    thin = thin,
    tune_burn = min(500L, nsim %/% 2L)
  )

  psc_fixed_beta <- as.numeric(as_draws(psc_fixed$draws)$beta_1)
  psc_full_beta <- as.numeric(as_draws(psc_full$draws)$beta_1)

  summ_row <- function(samples, label) {
    data.frame(
      source = label,
      n = length(samples),
      mean = mean(samples),
      sd = sd(samples),
      median = median(samples),
      hpd_lo = sim_hpd_interval(samples, 0.95)[["lower"]],
      hpd_hi = sim_hpd_interval(samples, 0.95)[["upper"]],
      stringsAsFactors = FALSE
    )
  }

  comparison <- rbind(
    summ_row(stan_beta, "Stan (fixed CFM)"),
    summ_row(psc_fixed_beta, "PSC fixed CFM (corrected MH)"),
    summ_row(psc_full_beta, "PSC full (CFM resampled)")
  )
  comparison <- transform(
    comparison,
    mean_bias = mean - SIM_DGP$beta,
    acc_rate = c(NA_real_, mean(psc_fixed$acc_trace, na.rm = TRUE), mean(psc_full$acc_trace, na.rm = TRUE)),
    proposal_scale = c(NA_real_, psc_fixed$proposal_scale_final, psc_full$proposal_scale_final)
  )

  cat("\nPSC MCMC validation (", overlap, ", ", sample_size, ", seed=", seed, ")\n", sep = "")
  cat("True log(HR) = ", round(SIM_DGP$beta, 4), "\n\n", sep = "")
  cmp_print <- comparison
  num_cols <- vapply(cmp_print, is.numeric, logical(1))
  cmp_print[num_cols] <- lapply(cmp_print[num_cols], function(x) round(x, 4))
  print(cmp_print, row.names = FALSE)

  ks_stan_psc <- stats::ks.test(stan_beta, psc_fixed_beta)
  cat("\nKS test (Stan vs PSC fixed CFM): D = ",
      round(ks_stan_psc$statistic, 4),
      ", p = ", format.pval(ks_stan_psc$p.value, digits = 3), "\n",
      sep = "")

  out <- list(
    overlap = overlap,
    sample_size = sample_size,
    seed = seed,
    true_beta = SIM_DGP$beta,
    comparison = comparison,
    stan_beta = stan_beta,
    psc_fixed_beta = psc_fixed_beta,
    psc_full_beta = psc_full_beta,
    ks_stan_psc = ks_stan_psc,
    lik_prep = lik_prep
  )

  if (save_results) {
    out_file <- file.path(paths$results_dir, "psc_mcmc_validation.rds")
    saveRDS(out, out_file)
    write.csv(comparison, file.path(paths$results_dir, "psc_mcmc_validation.csv"), row.names = FALSE)
    cat("\nSaved:", out_file, "\n")
  }

  invisible(out)
}

if (sys.nframe() == 0L || !interactive()) {
  PSC_SIM_DIR <- Sys.getenv("PSC_SIM_DIR", unset = normalizePath("."))
  assign("PSC_SIM_DIR", PSC_SIM_DIR, envir = .GlobalEnv)
  source(file.path(PSC_SIM_DIR, "setup.R"))
  validate_psc_mcmc()
}
