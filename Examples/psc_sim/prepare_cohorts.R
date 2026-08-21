# Shared cohort preparation for single-arm and hybrid analyses.

sim_prepare_single_arm <- function(contDat, trtDat) {
  contDat$time <- contDat$tmc
  contDat$status <- contDat$cenc
  trtDat$time <- trtDat$tm
  trtDat$status <- trtDat$cen
  contDat$trt <- 0
  trtDat$trt <- 1
  combDat <- rbind(contDat, trtDat)
  contDat$s.ob <- survival::Surv(contDat$time, contDat$status)
  trtDat$s.ob <- survival::Surv(trtDat$time, trtDat$status)
  combDat$s.ob <- survival::Surv(combDat$time, combDat$status)
  list(contDat = contDat, trtDat = trtDat, combDat = combDat)
}

sim_prepare_hybrid <- function(contDat, trtDat) {
  contDat$trt <- 0
  trtDat$trt <- rbinom(nrow(trtDat), 1, 0.5)
  contDat$time <- contDat$tmc
  contDat$status <- contDat$cenc
  trtDat$time <- trtDat$tm
  trtDat$status <- trtDat$cen
  trt_on_control <- trtDat$trt == 0
  trtDat$time[trt_on_control] <- trtDat$tmc[trt_on_control]
  trtDat$status[trt_on_control] <- trtDat$cenc[trt_on_control]
  contDat$source <- "hist"
  trtDat$source <- "cont"
  combDat <- rbind(contDat, trtDat)
  combDat$cen <- ifelse(combDat$source == "hist", combDat$cenc, combDat$cen)
  contDat$s.ob <- survival::Surv(contDat$time, contDat$status)
  trtDat$s.ob <- survival::Surv(trtDat$time, trtDat$status)
  combDat$s.ob <- survival::Surv(combDat$time, combDat$cen)
  list(contDat = contDat, trtDat = trtDat, combDat = combDat)
}

sim_factor_align <- function(contDat, trtDat) {
  contDat$X_3 <- factor(contDat$X_3)
  contDat$X_4 <- factor(contDat$X_4)
  contDat$X_5 <- factor(contDat$X_5)
  trtDat$X_3 <- factor(trtDat$X_3, levels = levels(contDat$X_3))
  trtDat$X_4 <- factor(trtDat$X_4, levels = levels(contDat$X_4))
  trtDat$X_5 <- factor(trtDat$X_5, levels = levels(contDat$X_5))
  trtDat$cen <- trtDat$status
  list(contDat = contDat, trtDat = trtDat)
}

sim_fit_cfm <- function(contDat) {
  cfm <- flexsurvspline(
    s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5,
    data = contDat,
    k = 3
  )
  pcfm <- flexsurvspline(
    s.ob ~ X_1 + X_3,
    data = contDat,
    k = 3
  )
  list(
    cfm = cfm,
    pcfm = pcfm,
    cfmw = as_pscCFM_dev(cfm),
    pcfmw = as_pscCFM_dev(pcfm)
  )
}
