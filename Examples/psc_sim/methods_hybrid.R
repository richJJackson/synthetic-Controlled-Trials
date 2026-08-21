# Hybrid trial methods for one simulation replicate.

sim_run_hybrid <- function(contDat, trtDat, defaults = SIM_DEFAULTS) {
  prep <- sim_prepare_hybrid(contDat, trtDat)
  contDat <- prep$contDat
  trtDat <- prep$trtDat
  combDat <- prep$combDat

  cm.rand <- coxph(s.ob ~ trt, data = trtDat)
  cm.rand.ad <- coxph(s.ob ~ X_1 + X_3 + trt, data = trtDat)
  cm.pooled <- coxph(s.ob ~ trt, data = combDat)
  cm.pooled.ad <- coxph(s.ob ~ X_1 + X_3 + trt, data = combDat)

  aligned <- sim_factor_align(contDat, trtDat)
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

  contCont <- combDat[which(combDat$trt == 0 & combDat$source == "cont"), ]
  synthDat <- combDat[-which(combDat$trt == 0 & combDat$source == "cont"), ]
  X <- as.matrix(synthDat[, c("X_1", "X_2", "X_3", "X_4", "X_5")])
  eb <- ebalance(Treatment = synthDat$trt, X = X)
  pX <- as.matrix(synthDat[, c("X_1", "X_3")])
  peb <- ebalance(Treatment = synthDat$trt, X = pX)
  synthDat$w <- 1
  synthDat$w[synthDat$trt == 0] <- eb$w
  synthDat$pw <- 1
  synthDat$pw[synthDat$trt == 0] <- peb$w
  contCont$w <- 1
  contCont$pw <- 1
  scData <- rbind(synthDat, contCont)
  scData$s.ob <- survival::Surv(scData$time, scData$cen)
  sc_cm <- coxph(s.ob ~ trt, data = scData, weights = scData$w)
  sc_pcm <- coxph(s.ob ~ trt, data = scData, weights = scData$pw)

  sr <- survreg(s.ob ~ 1, data = contDat, dist = "exponential")
  co <- as.numeric(sr$coefficients)
  se <- as.numeric(sqrt(sr$var))
  x <- as.numeric(trtDat$trt) - 1

  vague_baye <- sim_stan_fit(
    "exp",
    list(N = nrow(trtDat), time = trtDat$time, status = trtDat$status, x = x, lam0_mn = 0, lam0_t = 10),
    defaults
  )
  hybrid_baye <- sim_stan_fit(
    "exp",
    list(N = nrow(trtDat), time = trtDat$time, status = trtDat$status, x = x, lam0_mn = co, lam0_t = se),
    defaults
  )
  hybrid_cp_baye <- sim_stan_fit(
    "comm",
    list(N = nrow(trtDat), time = trtDat$time, status = trtDat$status, x = x, lam0_mn = co),
    defaults
  )

  CW_Data <- synthDat[which(synthDat$source == "hist"), ]
  hybrid_caseW_baye <- sim_stan_fit(
    "casew",
    list(
      N = nrow(trtDat), time = trtDat$time, status = trtDat$status, x = x,
      N0 = nrow(CW_Data), timec = CW_Data$time, statusc = CW_Data$cen, a0 = CW_Data$w
    ),
    defaults
  )
  hybrid_pcaseW_baye <- sim_stan_fit(
    "casew",
    list(
      N = nrow(trtDat), time = trtDat$time, status = trtDat$status, x = x,
      N0 = nrow(CW_Data), timec = CW_Data$time, statusc = CW_Data$cen, a0 = CW_Data$pw
    ),
    defaults
  )

  wald <- function(est, se) {
    z <- qnorm(0.975)
    list(lo = est - z * se, hi = est + z * se)
  }

  sim_replicate_out(
    est = c(
      summary(cm.rand)$coef[1, 1],
      summary(cm.rand.ad)$coef[3, 1],
      summary(cm.pooled)$coef[1, 1],
      summary(cm.pooled.ad)$coef[3, 1],
      psc_full_ci$est,
      psc_part_ci$est,
      summary(sc_cm)$coef[1, 1],
      summary(sc_pcm)$coef[1, 1],
      wpsc_full_ci$est,
      wpsc_part_ci$est,
      summary(vague_baye)$summary[2, 1],
      summary(hybrid_baye)$summary[2, 1],
      summary(hybrid_cp_baye)$summary[2, 1],
      summary(hybrid_caseW_baye)$summary[2, 1],
      summary(hybrid_pcaseW_baye)$summary[2, 1]
    ),
    se = c(
      summary(cm.rand)$coef[1, 3],
      summary(cm.rand.ad)$coef[3, 3],
      summary(cm.pooled)$coef[1, 3],
      summary(cm.pooled.ad)$coef[3, 3],
      psc_full_ci$se,
      psc_part_ci$se,
      summary(sc_cm)$coef[1, 3],
      summary(sc_pcm)$coef[1, 3],
      wpsc_full_ci$se,
      wpsc_part_ci$se,
      summary(vague_baye)$summary[2, 3],
      summary(hybrid_baye)$summary[2, 3],
      summary(hybrid_cp_baye)$summary[2, 3],
      summary(hybrid_caseW_baye)$summary[2, 3],
      summary(hybrid_pcaseW_baye)$summary[2, 3]
    ),
    lo = c(
      wald(summary(cm.rand)$coef[1, 1], summary(cm.rand)$coef[1, 3])$lo,
      wald(summary(cm.rand.ad)$coef[3, 1], summary(cm.rand.ad)$coef[3, 3])$lo,
      wald(summary(cm.pooled)$coef[1, 1], summary(cm.pooled)$coef[1, 3])$lo,
      wald(summary(cm.pooled.ad)$coef[3, 1], summary(cm.pooled.ad)$coef[3, 3])$lo,
      psc_full_ci$lo,
      psc_part_ci$lo,
      wald(summary(sc_cm)$coef[1, 1], summary(sc_cm)$coef[1, 3])$lo,
      wald(summary(sc_pcm)$coef[1, 1], summary(sc_pcm)$coef[1, 3])$lo,
      wpsc_full_ci$lo,
      wpsc_part_ci$lo,
      wald(summary(vague_baye)$summary[2, 1], summary(vague_baye)$summary[2, 3])$lo,
      wald(summary(hybrid_baye)$summary[2, 1], summary(hybrid_baye)$summary[2, 3])$lo,
      wald(summary(hybrid_cp_baye)$summary[2, 1], summary(hybrid_cp_baye)$summary[2, 3])$lo,
      wald(summary(hybrid_caseW_baye)$summary[2, 1], summary(hybrid_caseW_baye)$summary[2, 3])$lo,
      wald(summary(hybrid_pcaseW_baye)$summary[2, 1], summary(hybrid_pcaseW_baye)$summary[2, 3])$lo
    ),
    hi = c(
      wald(summary(cm.rand)$coef[1, 1], summary(cm.rand)$coef[1, 3])$hi,
      wald(summary(cm.rand.ad)$coef[3, 1], summary(cm.rand.ad)$coef[3, 3])$hi,
      wald(summary(cm.pooled)$coef[1, 1], summary(cm.pooled)$coef[1, 3])$hi,
      wald(summary(cm.pooled.ad)$coef[3, 1], summary(cm.pooled.ad)$coef[3, 3])$hi,
      psc_full_ci$hi,
      psc_part_ci$hi,
      wald(summary(sc_cm)$coef[1, 1], summary(sc_cm)$coef[1, 3])$hi,
      wald(summary(sc_pcm)$coef[1, 1], summary(sc_pcm)$coef[1, 3])$hi,
      wpsc_full_ci$hi,
      wpsc_part_ci$hi,
      wald(summary(vague_baye)$summary[2, 1], summary(vague_baye)$summary[2, 3])$hi,
      wald(summary(hybrid_baye)$summary[2, 1], summary(hybrid_baye)$summary[2, 3])$hi,
      wald(summary(hybrid_cp_baye)$summary[2, 1], summary(hybrid_cp_baye)$summary[2, 3])$hi,
      wald(summary(hybrid_caseW_baye)$summary[2, 1], summary(hybrid_caseW_baye)$summary[2, 3])$hi,
      wald(summary(hybrid_pcaseW_baye)$summary[2, 1], summary(hybrid_pcaseW_baye)$summary[2, 3])$hi
    )
  )
}
