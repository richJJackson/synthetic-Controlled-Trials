# Entropy-balance weights for historical controls.
# Weights sum to the trial sample size. Moments of X_1-X_6 in the weighted
# historical sample match the trial sample. Categorical covariates enter as
# level indicators, with level 1 as the reference.

entropy_balance_weights <- function(historical, trial) {
  both <- rbind(historical, trial)
  mm <- stats::model.matrix(
    ~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6,
    data = both
  )
  X <- mm[, -1, drop = FALSE]
  n0 <- nrow(historical)
  X0 <- X[seq_len(n0), , drop = FALSE]
  X1 <- X[-seq_len(n0), , drop = FALSE]
  target <- colMeans(X1)
  k <- ncol(X0)

  objective <- function(lambda) {
    u <- as.numeric(X0 %*% lambda)
    m <- max(u)
    logsum <- m + log(sum(exp(u - m)))
    as.numeric(logsum - sum(target * lambda))
  }
  gradient <- function(lambda) {
    u <- as.numeric(X0 %*% lambda)
    ew <- exp(u - max(u))
    w <- ew / sum(ew)
    as.numeric(colSums(X0 * w) - target)
  }

  fit <- stats::optim(rep(0, k), objective, gradient, method = "BFGS")
  u <- as.numeric(X0 %*% fit$par)
  ew <- pmax(exp(u - max(u)), 1e-12)
  w <- ew / sum(ew) * nrow(trial)
  imbalance <- max(abs(colSums(X0 * (w / sum(w))) - target))
  sum_w <- sum(w)
  list(
    w = w,
    converged = isTRUE(fit$convergence == 0) && imbalance < 0.05,
    max_imbalance = imbalance,
    sum_w = sum_w,
    ess = sum_w^2 / sum(w^2)
  )
}
