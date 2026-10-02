# Single-arm pooled Cox model with no covariate adjustment.
# Historical controls and treated trial patients are stacked.
# The estimate is a marginal log hazard ratio.

fit_pooled_unadjusted <- function(historical, trial) {
  historical$trt <- 0
  trial$trt <- 1
  comb <- rbind(historical, trial)
  comb$s.ob <- survival::Surv(comb$time, comb$status)
  fit <- survival::coxph(s.ob ~ trt, data = comb)
  co <- summary(fit)$coefficients
  est <- unname(co["trt", "coef"])
  se <- unname(co["trt", "se(coef)"])
  z <- stats::qnorm(0.975)
  # Every historical control enters with weight 1.
  data.frame(
    method = "Pooled unadjusted",
    est = est,
    se = se,
    lo = est - z * se,
    hi = est + z * se,
    max_imbalance = NA_real_,
    max_weight = 1,
    sum_w = nrow(historical),
    ess = nrow(historical),
    weight_on = "historical controls",
    accepted = NA,
    calibration_slope = NA_real_,
    stringsAsFactors = FALSE
  )
}
