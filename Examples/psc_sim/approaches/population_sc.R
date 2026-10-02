# Population synthetic control.
# Entropy-balance the historical controls to the trial covariate distribution,
# then fit a weighted Cox model of treatment only.

fit_population_sc <- function(historical, trial) {
  bal <- entropy_balance_weights(historical, trial)
  historical$trt <- 0
  trial$trt <- 1
  historical$w <- bal$w
  trial$w <- 1
  comb <- rbind(historical, trial)
  comb$s.ob <- survival::Surv(comb$time, comb$status)
  fit <- survival::coxph(s.ob ~ trt, data = comb, weights = w)
  co <- summary(fit)$coefficients
  est <- unname(co["trt", "coef"])
  se <- unname(co["trt", "se(coef)"])
  z <- stats::qnorm(0.975)
  data.frame(
    method = "Population synthetic control",
    est = est,
    se = se,
    lo = est - z * se,
    hi = est + z * se,
    max_imbalance = bal$max_imbalance,
    max_weight = max(bal$w),
    sum_w = bal$sum_w,
    ess = bal$ess,
    weight_on = "historical controls",
    accepted = NA,
    calibration_slope = NA_real_,
    stringsAsFactors = FALSE
  )
}
