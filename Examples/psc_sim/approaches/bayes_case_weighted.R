# Case-weighted Bayesian historical control for a single-arm trial.
# The historical likelihood for the control log-mean is weighted by the
# entropy-balance weights. The trial likelihood is the treated hazard
# exp(-lambda0 + beta). Priors are N(0, 10) on both parameters.
# Posterior summaries use a random-walk Metropolis sampler.

fit_bayes_case_weighted <- function(historical, trial, nsim = 8000L, burn = 2000L, seed = 1L) {
  bal <- entropy_balance_weights(historical, trial)
  a0 <- bal$w
  log_post <- function(theta) {
    lambda0 <- theta[1]
    beta <- theta[2]
    lp <- stats::dnorm(lambda0, 0, 10, log = TRUE) +
      stats::dnorm(beta, 0, 10, log = TRUE)
    rate_t <- exp(-lambda0 + beta)
    rate_c <- exp(-lambda0)
    lp <- lp + sum(trial$status * log(rate_t) - rate_t * trial$time)
    lp <- lp + sum(a0 * (historical$status * log(rate_c) - rate_c * historical$time))
    lp
  }

  set.seed(seed)
  theta <- c(0, 0)
  cur <- log_post(theta)
  step <- c(0.15, 0.25)
  draws <- matrix(NA_real_, nsim, 2)
  accepts <- 0L
  for (i in seq_len(nsim)) {
    prop <- theta + stats::rnorm(2, 0, step)
    pr <- log_post(prop)
    if (is.finite(pr) && log(stats::runif(1)) < (pr - cur)) {
      theta <- prop
      cur <- pr
      accepts <- accepts + 1L
    }
    draws[i, ] <- theta
    if (i %% 500L == 0L && i < burn) {
      rate <- accepts / 500
      if (rate < 0.15) step <- step * 0.7
      if (rate > 0.45) step <- step * 1.3
      accepts <- 0L
    }
  }
  beta <- draws[(burn + 1L):nsim, 2]
  q <- stats::quantile(beta, c(0.025, 0.975), names = FALSE)
  data.frame(
    method = "Bayesian case-weighted",
    est = stats::median(beta),
    se = stats::sd(beta),
    lo = q[1],
    hi = q[2],
    max_imbalance = bal$max_imbalance,
    max_weight = max(a0),
    sum_w = bal$sum_w,
    ess = bal$ess,
    weight_on = "historical controls",
    accepted = NA,
    calibration_slope = NA_real_,
    stringsAsFactors = FALSE
  )
}
