# Patient-level LEAP exchange probabilities against covariate closeness to the trial.
# Same samples as hybrid_100.R. LEAP only.

sim_dir <- normalizePath("~/Documents/GitHub/synthetic-Controlled-Trials/Examples/psc_sim", mustWork = TRUE)
load(file.path(sim_dir, "data/master.RData"))
done <- utils::read.csv(file.path(sim_dir, "data/hybrid_100.csv"))

overlaps <- c(dataset_1 = "Same margins", dataset_2 = "Moderate overlap", dataset_3 = "Poor overlap")
seeds <- c(dataset_1 = 11L, dataset_2 = 12L, dataset_3 = 13L)

picks <- do.call(rbind, lapply(names(overlaps), function(nm) {
  sub <- done[done$overlap == overlaps[[nm]], ]
  med_target <- stats::median(sub$n_exchangeable)
  data.frame(
    nm = nm,
    overlap = overlaps[[nm]],
    role = c("few_exchangeable", "typical", "most_exchangeable"),
    rep = c(
      sub$rep[which.min(sub$n_exchangeable)],
      sub$rep[which.min(abs(sub$n_exchangeable - med_target))],
      sub$rep[which.max(sub$n_exchangeable)]
    ),
    stringsAsFactors = FALSE
  )
}))

arm_of <- function(dat, outcome) {
  if (outcome == "control") {
    dat$time <- dat$t0
    dat$status <- dat$event0
  } else {
    dat$time <- dat$t1
    dat$status <- dat$event1
  }
  dat
}

patient_rows <- list()
example_rows <- list()

for (k in seq_len(nrow(picks))) {
  nm <- picks$nm[k]
  rep_i <- picks$rep[k]
  seed_h <- unname(seeds[[nm]]) + (rep_i - 1L) * 1000L
  set.seed(seed_h)
  historical <- arm_of(master$dataset_1[sample.int(nrow(master$dataset_1), 500L), ], "control")
  set.seed(seed_h + 100L)
  trial_source <- master[[nm]]
  trial <- trial_source[sample.int(nrow(trial_source), 75L), ]
  set.seed(seed_h + 200L)
  trial$arm <- sample(c(rep(0L, 38L), rep(1L, 37L)))
  trial$time <- ifelse(trial$arm == 0L, trial$t0, trial$t1)
  trial$status <- ifelse(trial$arm == 0L, trial$event0, trial$event1)
  for (v in c("X_4", "X_5", "X_6")) {
    historical[[v]] <- factor(historical[[v]], levels = 1:3)
    trial[[v]] <- factor(trial[[v]], levels = 1:3)
  }

  sr <- survival::survreg(
    survival::Surv(time, status) ~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6,
    data = historical,
    dist = "exponential"
  )
  mu_hist <- unname(stats::coef(sr)["(Intercept)"])
  X_trial <- stats::model.matrix(~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6, data = trial)[, -1, drop = FALSE]
  X_hist <- stats::model.matrix(~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6, data = historical)[, -1, drop = FALSE]

  set.seed(seed_h + 4L)
  n_h <- nrow(historical)
  gamma1 <- 0.8
  lambda2 <- mu_hist
  theta_l <- c(mu_hist, -unname(stats::coef(sr)[-1]), log(0.7))
  loglik_class <- function(lambda0, gamma, lambda2) {
    eta1 <- -lambda0 + as.numeric(X_hist %*% gamma)
    ll1 <- historical$status * eta1 - exp(eta1) * historical$time
    eta2 <- -lambda2
    ll2 <- historical$status * eta2 - exp(eta2) * historical$time
    cbind(ll1, ll2)
  }
  log_post_leap <- function(theta, class1, lambda2) {
    lambda0 <- theta[1]
    gamma <- theta[2:10]
    beta <- theta[11]
    lp <- stats::dnorm(lambda0, 0, 10, log = TRUE)
    lp <- lp + sum(stats::dnorm(gamma, 0, 10, log = TRUE))
    lp <- lp + stats::dnorm(beta, 0, 10, log = TRUE)
    lp <- lp + stats::dnorm(lambda2, 0, 10, log = TRUE)
    eta_t <- -lambda0 + as.numeric(X_trial %*% gamma) + beta * trial$arm
    lp <- lp + sum(trial$status * eta_t - exp(eta_t) * trial$time)
    if (any(class1)) {
      eta1 <- -lambda0 + as.numeric(X_hist[class1, , drop = FALSE] %*% gamma)
      lp <- lp + sum(historical$status[class1] * eta1 - exp(eta1) * historical$time[class1])
    }
    lp
  }
  class1 <- rep(TRUE, n_h)
  cur <- log_post_leap(theta_l, class1, lambda2)
  step_l <- rep(0.06, 11L)
  step2 <- 0.15
  n_leap <- 4000L
  burn_l <- 1000L
  p_exch <- matrix(NA_real_, n_leap, n_h)
  lambda2_draw <- rep(NA_real_, n_leap)
  theta_draw <- matrix(NA_real_, n_leap, 11L)
  acc_l <- 0L
  for (i in seq_len(n_leap)) {
    ll <- loglik_class(theta_l[1], theta_l[2:10], lambda2)
    logp1 <- log(gamma1) + ll[, 1]
    logp2 <- log(1 - gamma1) + ll[, 2]
    m <- pmax(logp1, logp2)
    p1 <- exp(logp1 - m) / (exp(logp1 - m) + exp(logp2 - m))
    class1 <- stats::runif(n_h) < p1
    p_exch[i, ] <- p1
    n1 <- sum(class1)
    gamma1 <- stats::rbeta(1, 1 + n1, 1 + (n_h - n1))
    cur <- log_post_leap(theta_l, class1, lambda2)
    prop <- theta_l + stats::rnorm(11L, 0, step_l)
    pr <- log_post_leap(prop, class1, lambda2)
    if (is.finite(pr) && log(stats::runif(1)) < (pr - cur)) {
      theta_l <- prop
      cur <- pr
      acc_l <- acc_l + 1L
    }
    if (any(!class1)) {
      prop2 <- lambda2 + stats::rnorm(1, 0, step2)
      eta_old <- -lambda2
      eta_new <- -prop2
      use <- !class1
      ll_old <- sum(historical$status[use] * eta_old - exp(eta_old) * historical$time[use]) +
        stats::dnorm(lambda2, 0, 10, log = TRUE)
      ll_new <- sum(historical$status[use] * eta_new - exp(eta_new) * historical$time[use]) +
        stats::dnorm(prop2, 0, 10, log = TRUE)
      if (is.finite(ll_new) && log(stats::runif(1)) < (ll_new - ll_old)) {
        lambda2 <- prop2
        cur <- log_post_leap(theta_l, class1, lambda2)
      }
    }
    theta_draw[i, ] <- theta_l
    lambda2_draw[i] <- lambda2
    if (i %% 200L == 0L && i < burn_l) {
      rate <- acc_l / 200
      if (rate < 0.15) step_l <- step_l * 0.75
      if (rate > 0.35) step_l <- step_l * 1.1
      acc_l <- 0L
    }
  }
  keep <- (burn_l + 1L):n_leap
  p_bar <- colMeans(p_exch[keep, , drop = FALSE])
  theta_bar <- colMeans(theta_draw[keep, , drop = FALSE])
  lambda2_bar <- mean(lambda2_draw[keep])
  ll_bar <- loglik_class(theta_bar[1], theta_bar[2:10], lambda2_bar)
  gap <- ll_bar[, 1] - ll_bar[, 2]

  eig <- eigen(stats::cov(X_trial), symmetric = TRUE)
  eig$values[eig$values < 1e-8] <- 1e-8
  S <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
  center <- colMeans(X_trial)
  d_h <- stats::mahalanobis(X_hist, center, S)
  d_t <- stats::mahalanobis(X_trial, center, S)
  closer_than_trial <- vapply(d_h, function(d) mean(d_t >= d), numeric(1))

  pts <- data.frame(
    overlap = picks$overlap[k],
    role = picks$role[k],
    rep = rep_i,
    id = historical$id,
    p_exchangeable = p_bar,
    closer_than_trial = closer_than_trial,
    ll_gap = gap,
    time = historical$time,
    event = historical$status,
    X_1 = historical$X_1,
    X_2 = historical$X_2,
    X_3 = historical$X_3,
    X_4 = as.integer(as.character(historical$X_4)),
    X_5 = as.integer(as.character(historical$X_5)),
    X_6 = as.integer(as.character(historical$X_6)),
    stringsAsFactors = FALSE
  )
  patient_rows[[k]] <- pts
  good <- pts[pts$closer_than_trial >= 0.5, ]
  good <- good[order(good$p_exchangeable), ]
  take <- head(good, 4)
  example_rows[[k]] <- take
  cat(
    picks$overlap[k], picks$role[k], "rep", rep_i,
    "median p", round(stats::median(p_bar), 3),
    "p<0.5", sum(p_bar < 0.5),
    "good match", nrow(good),
    "good match and p<0.5", sum(good$p_exchangeable < 0.5),
    "cor(p, closeness)", round(stats::cor(p_bar, closer_than_trial), 3),
    "\n"
  )
  flush.console()
}

patients <- do.call(rbind, patient_rows)
examples <- do.call(rbind, example_rows)
utils::write.csv(patients, file.path(sim_dir, "data/leap_match_patients.csv"), row.names = FALSE)
utils::write.csv(examples, file.path(sim_dir, "data/leap_match_examples.csv"), row.names = FALSE)
cat("Saved patient file", nrow(patients), "example rows", nrow(examples), "\n")
print(examples[, c("overlap", "role", "rep", "p_exchangeable", "closer_than_trial", "ll_gap", "time", "event", "X_1", "X_2", "X_3", "X_4", "X_5", "X_6")])
