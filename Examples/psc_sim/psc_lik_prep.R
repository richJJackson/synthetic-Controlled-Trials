# Precompute Royston-Parmar likelihood pieces with fixed CFM coefficients.

psc_lik_prep_fixed <- function(pscOb, co = NULL) {
  co <- co %||% pscOb$co
  haz_co <- co[names(co) %in% names(pscOb$haz_co)]
  cov_co <- co[names(co) %in% names(pscOb$cov_co)]

  lam <- pscOb$lam
  kn <- pscOb$kn
  k <- pscOb$k

  time <- pscOb$DC$Y$time
  cen <- pscOb$DC$Y$cen
  cov <- pscOb$DC$X

  logt <- log(time)
  lp <- as.numeric(cov %*% cov_co)

  z <- NULL
  z_h <- NULL
  for (i in seq_len(k)) {
    zt <- modp(logt - kn[(i + 1)])^3 -
      lam[(i + 1)] * modp(logt - kn[1])^3 -
      (1 - lam[(i + 1)]) * modp(logt - kn[length(kn)])^3
    z <- cbind(z, zt)
    zt_h <- (modp(logt - kn[(i + 1)])^2 -
      lam[(i + 1)] * modp(logt - kn[1])^2 -
      (1 - lam[(i + 1)]) * modp(logt - kn[length(kn)])^2)
    z_h <- cbind(z_h, zt_h)
  }

  H0 <- exp(haz_co[1] + haz_co[2] * logt + z %*% haz_co[3:(2 + k)])
  h0 <- (H0 / time) * (haz_co[2] + 3 * z_h %*% haz_co[3:(2 + k)])

  data.frame(
    H0 = as.numeric(H0),
    h0 = as.numeric(h0),
    lp = lp,
    status = as.integer(cen),
    stringsAsFactors = FALSE
  )
}

psc_loglik_fixed_beta <- function(beta, prep) {
  H <- prep$H0 * exp(prep$lp + beta)
  h <- prep$h0 * exp(prep$lp + beta)
  S <- exp(-H)
  ll <- prep$status * log(S * h + 1e-16) + (1 - prep$status) * log(S + 1e-16)
  sum(ll)
}
