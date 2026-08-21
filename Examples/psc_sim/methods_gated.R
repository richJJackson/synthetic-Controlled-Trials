# Single-arm replicate: CFM internal-validation gate + optional PSC (doubled MCMC).

GATED_SA_METHOD_NAMES <- c(
  "Pooled (full adj)", "Pooled (partial adj)",
  "PSC full", "PSC partial",
  "PSC full (pw)", "PSC partial (pw)"
)

#' One gated single-arm replicate.
#'
#' Fits CFMs, runs internal validation on controls, and only runs PSC if the
#' gate passes. PSC uses defaults$psc_nsim (expected 10000 for this experiment).
#' Always returns pooled Cox references for comparison.
sim_run_single_arm_gated <- function(
    contDat,
    trtDat,
    defaults = SIM_DEFAULTS,
    cal_tol = 0.05,
    cindex_min = 0.55,
    validate_seed = 1L) {

  prep <- sim_prepare_single_arm(contDat, trtDat)
  contDat <- prep$contDat
  trtDat <- prep$trtDat
  combDat <- prep$combDat

  cm_adj <- coxph(s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5 + trt, data = combDat)
  cm_padj <- coxph(s.ob ~ X_1 + X_3 + trt, data = combDat)
  pool_full <- sim_wald_ci(coef(cm_adj)[["trt"]], sqrt(vcov(cm_adj)["trt", "trt"]))
  pool_part <- sim_wald_ci(coef(cm_padj)[["trt"]], sqrt(vcov(cm_padj)["trt", "trt"]))

  aligned <- sim_factor_align(contDat, trtDat)
  contDat <- aligned$contDat
  trtDat <- aligned$trtDat

  gate_full <- sim_cfm_internal_validate(
    contDat,
    formula = s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5,
    cal_tol = cal_tol,
    cindex_min = cindex_min,
    seed = validate_seed
  )
  gate_part <- sim_cfm_internal_validate(
    contDat,
    formula = s.ob ~ X_1 + X_3,
    cal_tol = cal_tol,
    cindex_min = cindex_min,
    seed = validate_seed + 1L
  )

  # Gate on the full CFM for PSC full / pw-full; partial CFM for PSC partial.
  na_ci <- list(est = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_)
  psc_full_ci <- na_ci
  psc_part_ci <- na_ci
  wpsc_full_ci <- na_ci
  wpsc_part_ci <- na_ci

  models <- NULL
  if (gate_full$pass || gate_part$pass) {
    models <- sim_fit_cfm(contDat)
  }

  if (gate_full$pass) {
    pw <- proximityWeights(models$cfmw, trtDat, method = "joint")
    psc_full <- pscfit(
      models$cfmw, trtDat,
      nsim = defaults$psc_nsim, nchain = defaults$psc_nchain
    )
    wpsc_full <- pscfit(
      models$cfmw, trtDat, proximity_weights = pw,
      nsim = defaults$psc_nsim, nchain = defaults$psc_nchain
    )
    psc_full_ci <- sim_psc_beta_interval(psc_full)
    wpsc_full_ci <- sim_psc_beta_interval(wpsc_full)
  }

  if (gate_part$pass) {
    ppw <- proximityWeights(models$pcfmw, trtDat, method = "joint")
    psc_part <- pscfit(
      models$pcfmw, trtDat,
      nsim = defaults$psc_nsim, nchain = defaults$psc_nchain
    )
    wpsc_part <- pscfit(
      models$pcfmw, trtDat, proximity_weights = ppw,
      nsim = defaults$psc_nsim, nchain = defaults$psc_nchain
    )
    psc_part_ci <- sim_psc_beta_interval(psc_part)
    wpsc_part_ci <- sim_psc_beta_interval(wpsc_part)
  }

  out <- sim_replicate_out(
    est = c(
      pool_full$est, pool_part$est,
      psc_full_ci$est, psc_part_ci$est,
      wpsc_full_ci$est, wpsc_part_ci$est
    ),
    se = c(
      pool_full$se, pool_part$se,
      psc_full_ci$se, psc_part_ci$se,
      wpsc_full_ci$se, wpsc_part_ci$se
    ),
    lo = c(
      pool_full$lo, pool_part$lo,
      psc_full_ci$lo, psc_part_ci$lo,
      wpsc_full_ci$lo, wpsc_part_ci$lo
    ),
    hi = c(
      pool_full$hi, pool_part$hi,
      psc_full_ci$hi, psc_part_ci$hi,
      wpsc_full_ci$hi, wpsc_part_ci$hi
    )
  )

  list(
    res = out,
    gate_full_pass = gate_full$pass,
    gate_part_pass = gate_part$pass,
    gate_full_reason = gate_full$reason,
    gate_part_reason = gate_part$reason,
    gate_full_cal_abs = gate_full$cal_abs,
    gate_part_cal_abs = gate_part$cal_abs,
    gate_full_cindex = gate_full$cindex,
    gate_part_cindex = gate_part$cindex
  )
}
