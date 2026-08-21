# Single-arm trial methods for one simulation replicate.

sim_run_single_arm <- function(contDat, trtDat, defaults = SIM_DEFAULTS) {
  prep <- sim_prepare_single_arm(contDat, trtDat)
  contDat <- prep$contDat
  trtDat <- prep$trtDat
  combDat <- prep$combDat

  cm <- coxph(s.ob ~ trt, data = combDat)
  cm_adj <- coxph(s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5 + trt, data = combDat)
  cm_padj <- coxph(s.ob ~ X_1 + X_3 + trt, data = combDat)

  aligned <- sim_factor_align(contDat, trtDat)
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

  X <- as.matrix(combDat[, c("X_1", "X_2", "X_3", "X_4")])
  eb <- ebalance(Treatment = combDat$trt, X = X)
  combDat$w <- 1
  combDat$w[combDat$trt == 0] <- eb$w
  pX <- as.matrix(combDat[, c("X_1", "X_3")])
  peb <- ebalance(Treatment = combDat$trt, X = pX)
  combDat$pw <- 1
  combDat$pw[combDat$trt == 0] <- peb$w
  sc_cm <- coxph(s.ob ~ trt, data = combDat, weights = combDat$w)
  sc_pcm <- coxph(s.ob ~ trt, data = combDat, weights = combDat$pw)

  sr <- survreg(s.ob ~ 1, data = contDat, dist = "exponential")
  co <- as.numeric(sr$coefficients)
  se <- as.numeric(sqrt(sr$var))
  singArm_baye <- sim_stan_fit(
    "exp",
    list(
      N = nrow(trtDat),
      time = trtDat$time,
      status = trtDat$status,
      x = trtDat$trt,
      lam0_mn = co,
      lam0_t = se
    ),
    defaults
  )
  cp_w <- combDat$w[combDat$trt == 0]
  singArm_caseW_baye <- sim_stan_fit(
    "casew",
    list(
      N = nrow(trtDat),
      time = trtDat$time,
      status = trtDat$status,
      x = trtDat$trt,
      N0 = nrow(contDat),
      timec = contDat$time,
      statusc = contDat$status,
      a0 = cp_w
    ),
    defaults
  )

  sim_replicate_out(
    est = c(
      summary(cm)$coef[1, 1],
      summary(cm_adj)$coef[6, 1],
      summary(cm_padj)$coef[3, 1],
      psc_full_ci$est,
      psc_part_ci$est,
      summary(sc_cm)$coef[1, 1],
      summary(sc_pcm)$coef[1, 1],
      summary(singArm_baye)$summary[2, 1],
      summary(singArm_caseW_baye)$summary[2, 1],
      wpsc_full_ci$est,
      wpsc_part_ci$est
    ),
    se = c(
      summary(cm)$coef[1, 3],
      summary(cm_adj)$coef[6, 3],
      summary(cm_padj)$coef[3, 3],
      psc_full_ci$se,
      psc_part_ci$se,
      summary(sc_cm)$coef[1, 3],
      summary(sc_pcm)$coef[1, 3],
      summary(singArm_baye)$summary[2, 3],
      summary(singArm_caseW_baye)$summary[2, 3],
      wpsc_full_ci$se,
      wpsc_part_ci$se
    ),
    lo = c(
      summary(cm)$coef[1, 1] - qnorm(0.975) * summary(cm)$coef[1, 3],
      summary(cm_adj)$coef[6, 1] - qnorm(0.975) * summary(cm_adj)$coef[6, 3],
      summary(cm_padj)$coef[3, 1] - qnorm(0.975) * summary(cm_padj)$coef[3, 3],
      psc_full_ci$lo,
      psc_part_ci$lo,
      summary(sc_cm)$coef[1, 1] - qnorm(0.975) * summary(sc_cm)$coef[1, 3],
      summary(sc_pcm)$coef[1, 1] - qnorm(0.975) * summary(sc_pcm)$coef[1, 3],
      summary(singArm_baye)$summary[2, 1] - qnorm(0.975) * summary(singArm_baye)$summary[2, 3],
      summary(singArm_caseW_baye)$summary[2, 1] - qnorm(0.975) * summary(singArm_caseW_baye)$summary[2, 3],
      wpsc_full_ci$lo,
      wpsc_part_ci$lo
    ),
    hi = c(
      summary(cm)$coef[1, 1] + qnorm(0.975) * summary(cm)$coef[1, 3],
      summary(cm_adj)$coef[6, 1] + qnorm(0.975) * summary(cm_adj)$coef[6, 3],
      summary(cm_padj)$coef[3, 1] + qnorm(0.975) * summary(cm_padj)$coef[3, 3],
      psc_full_ci$hi,
      psc_part_ci$hi,
      summary(sc_cm)$coef[1, 1] + qnorm(0.975) * summary(sc_cm)$coef[1, 3],
      summary(sc_pcm)$coef[1, 1] + qnorm(0.975) * summary(sc_pcm)$coef[1, 3],
      summary(singArm_baye)$summary[2, 1] + qnorm(0.975) * summary(singArm_baye)$summary[2, 3],
      summary(singArm_caseW_baye)$summary[2, 1] + qnorm(0.975) * summary(singArm_caseW_baye)$summary[2, 3],
      wpsc_full_ci$hi,
      wpsc_part_ci$hi
    )
  )
}
