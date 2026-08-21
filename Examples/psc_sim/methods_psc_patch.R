# PSC-only replicate runners for patching existing simulation results.

sim_psc_method_rows <- function(trial_type = c("single_arm", "hybrid")) {
  trial_type <- match.arg(trial_type)
  method_names <- if (trial_type == "single_arm") SA_METHOD_NAMES else HYBRID_METHOD_NAMES
  psc_names <- if (trial_type == "single_arm") {
    c("PSC full", "PSC partial", "PSC full (pw)", "PSC partial (pw)")
  } else {
    c(
      "PSC full (combined)", "PSC partial (combined)",
      "PSC full (pw, combined)", "PSC partial (pw, combined)"
    )
  }
  idx <- match(psc_names, method_names)
  if (any(is.na(idx))) {
    stop("PSC method names not found in method list.")
  }
  idx
}

sim_run_psc_single_arm <- function(contDat, trtDat, defaults = SIM_DEFAULTS) {
  prep <- sim_prepare_single_arm(contDat, trtDat)
  aligned <- sim_factor_align(prep$contDat, prep$trtDat)
  contDat <- aligned$contDat
  trtDat <- aligned$trtDat
  models <- sim_fit_cfm(contDat)

  pw <- proximityWeights(models$cfmw, trtDat, method = "joint")
  ppw <- proximityWeights(models$pcfmw, trtDat, method = "joint")

  psc_full <- pscfit(
    models$cfmw, trtDat, nsim = defaults$psc_nsim, nchain = defaults$psc_nchain
  )
  psc_part <- pscfit(
    models$pcfmw, trtDat, nsim = defaults$psc_nsim, nchain = defaults$psc_nchain
  )
  wpsc_full <- pscfit(
    models$cfmw, trtDat, proximity_weights = pw,
    nsim = defaults$psc_nsim, nchain = defaults$psc_nchain
  )
  wpsc_part <- pscfit(
    models$pcfmw, trtDat, proximity_weights = ppw,
    nsim = defaults$psc_nsim, nchain = defaults$psc_nchain
  )

  psc_full_ci <- sim_psc_beta_interval(psc_full)
  psc_part_ci <- sim_psc_beta_interval(psc_part)
  wpsc_full_ci <- sim_psc_beta_interval(wpsc_full)
  wpsc_part_ci <- sim_psc_beta_interval(wpsc_part)

  sim_replicate_out(
    est = c(psc_full_ci$est, psc_part_ci$est, wpsc_full_ci$est, wpsc_part_ci$est),
    se = c(psc_full_ci$se, psc_part_ci$se, wpsc_full_ci$se, wpsc_part_ci$se),
    lo = c(psc_full_ci$lo, psc_part_ci$lo, wpsc_full_ci$lo, wpsc_part_ci$lo),
    hi = c(psc_full_ci$hi, psc_part_ci$hi, wpsc_full_ci$hi, wpsc_part_ci$hi)
  )
}

sim_run_psc_hybrid <- function(contDat, trtDat, defaults = SIM_DEFAULTS) {
  prep <- sim_prepare_hybrid(contDat, trtDat)
  aligned <- sim_factor_align(prep$contDat, prep$trtDat)
  contDat <- aligned$contDat
  trtDat <- aligned$trtDat
  trtDat$trt <- factor(trtDat$trt + 1)
  models <- sim_fit_cfm(contDat)

  pw <- proximityWeights(models$cfmw, trtDat, method = "joint")
  ppw <- proximityWeights(models$pcfmw, trtDat, method = "joint")

  psc_full <- pscfit(
    models$cfmw, trtDat, nsim = defaults$psc_nsim,
    nchain = defaults$psc_nchain, trt = trtDat$trt
  )
  psc_part <- pscfit(
    models$pcfmw, trtDat, nsim = defaults$psc_nsim,
    nchain = defaults$psc_nchain, trt = trtDat$trt
  )
  wpsc_full <- pscfit(
    models$cfmw, trtDat, proximity_weights = pw,
    nsim = defaults$psc_nsim, nchain = defaults$psc_nchain, trt = trtDat$trt
  )
  wpsc_part <- pscfit(
    models$pcfmw, trtDat, proximity_weights = ppw,
    nsim = defaults$psc_nsim, nchain = defaults$psc_nchain, trt = trtDat$trt
  )

  psc_full_ci <- sim_psc_sample_interval(as.numeric(unlist(pscComb(psc_full))))
  psc_part_ci <- sim_psc_sample_interval(as.numeric(unlist(pscComb(psc_part))))
  wpsc_full_ci <- sim_psc_sample_interval(as.numeric(unlist(pscComb(wpsc_full))))
  wpsc_part_ci <- sim_psc_sample_interval(as.numeric(unlist(pscComb(wpsc_part))))

  sim_replicate_out(
    est = c(psc_full_ci$est, psc_part_ci$est, wpsc_full_ci$est, wpsc_part_ci$est),
    se = c(psc_full_ci$se, psc_part_ci$se, wpsc_full_ci$se, wpsc_part_ci$se),
    lo = c(psc_full_ci$lo, psc_part_ci$lo, wpsc_full_ci$lo, wpsc_part_ci$lo),
    hi = c(psc_full_ci$hi, psc_part_ci$hi, wpsc_full_ci$hi, wpsc_part_ci$hi)
  )
}
