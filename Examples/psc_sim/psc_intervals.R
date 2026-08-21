# Posterior intervals for PSC methods in the simulation study.

#' Shortest-interval 95% HPD for a numeric sample (Chen-Shao algorithm).
sim_hpd_interval <- function(x, prob = 0.95) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  if (n < 2L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  k <- ceiling(prob * n)
  if (k >= n) {
    k <- n - 1L
  }
  if (k < 1L) {
    k <- 1L
  }
  span <- x[(1L + k):n] - x[1L:(n - k)]
  best <- which.min(span)
  c(lower = x[best], upper = x[best + k])
}

sim_psc_beta_interval <- function(psc_fit, prob = 0.95) {
  drs <- psc_fit$draws
  if (is.null(drs)) {
    stop("PSC fit has no posterior draws.")
  }
  beta_var <- "beta_1"
  if (!beta_var %in% names(drs)) {
    beta_vars <- grep("^beta_", names(drs), value = TRUE)
    if (!length(beta_vars)) {
      stop("No beta draws found in PSC fit.")
    }
    beta_var <- beta_vars[1]
  }
  samples <- as.numeric(drs[[beta_var]])
  sim_psc_sample_interval(samples, prob = prob)
}

sim_psc_sample_interval <- function(samples, prob = 0.95) {
  samples <- as.numeric(samples)
  samples <- samples[is.finite(samples)]
  if (length(samples) < 2L) {
    return(list(est = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  hpd <- sim_hpd_interval(samples, prob = prob)
  list(
    est = stats::median(samples),
    se = stats::sd(samples),
    lo = unname(hpd[["lower"]]),
    hi = unname(hpd[["upper"]])
  )
}

sim_wald_ci <- function(est, se, z = qnorm(0.975)) {
  list(
    est = est,
    se = se,
    lo = est - z * se,
    hi = est + z * se
  )
}

sim_replicate_out <- function(est, se, lo = NULL, hi = NULL, z = qnorm(0.975)) {
  if (is.null(lo) || is.null(hi)) {
    lo <- est - z * se
    hi <- est + z * se
  }
  cbind(est = est, se = se, lo = lo, hi = hi)
}

sim_expand_res_array <- function(res.array, n_methods = NULL, nsim = NULL) {
  if (is.null(res.array) || length(res.array) == 0L) {
    return(res.array)
  }
  n_col <- dim(res.array)[2]
  if (n_col >= SIM_RES_NCOL) {
    return(res.array)
  }
  if (is.null(n_methods)) {
    n_methods <- dim(res.array)[1]
  }
  if (is.null(nsim)) {
    nsim <- dim(res.array)[3]
  }
  out <- array(NA_real_, dim = c(n_methods, SIM_RES_NCOL, nsim))
  out[, seq_len(n_col), seq_len(nsim)] <- res.array
  if (n_col >= 2L) {
    z <- qnorm(0.975)
    est <- res.array[, 1, , drop = FALSE]
    se <- res.array[, 2, , drop = FALSE]
    out[, 3, ] <- est[, 1, ] - z * se[, 1, ]
    out[, 4, ] <- est[, 1, ] + z * se[, 1, ]
  }
  out
}
