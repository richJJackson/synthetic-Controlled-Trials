# Hybrid simulations, same margins. Replicates 11 to 20, appended to the first ten.
# 500 historical controls from dataset 1.
# 75 trial patients from dataset 1, randomized 38:37.
# Control outcomes are t0. Treated outcomes are t1.

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

sim_dir <- normalizePath("~/Documents/GitHub/synthetic-Controlled-Trials/Examples/psc_sim", mustWork = TRUE)
source(file.path(sim_dir, "approaches/entropy_balance.R"))
source(file.path(sim_dir, "approaches/psc.R"))
load(file.path(sim_dir, "data/master.RData"))

load_psc_for_prototype()
assignInNamespace("txtProgressBar", function(...) structure(list(), class = "txtProgressBar"), "utils")
assignInNamespace("setTxtProgressBar", function(pb, ...) invisible(pb), "utils")

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

run_rep <- function(seed_h) {
set.seed(seed_h)
historical <- arm_of(master$dataset_1[sample.int(nrow(master$dataset_1), 500L), ], "control")
set.seed(seed_h + 100L)
trial_rows <- master$dataset_1[sample.int(nrow(master$dataset_1), 75L), ]
set.seed(seed_h + 200L)
trial_rows$arm <- sample(c(rep(0L, 38L), rep(1L, 37L)))
trial <- trial_rows
trial$time <- ifelse(trial$arm == 0L, trial$t0, trial$t1)
trial$status <- ifelse(trial$arm == 0L, trial$event0, trial$event1)
trial$trt <- trial$arm

for (nm in c("X_4", "X_5", "X_6")) {
  historical[[nm]] <- factor(historical[[nm]], levels = 1:3)
  trial[[nm]] <- factor(trial[[nm]], levels = 1:3)
}
historical$s.ob <- survival::Surv(historical$time, historical$status)
trial$s.ob <- survival::Surv(trial$time, trial$status)
trial$cen <- trial$status

cohort_line <- function(dat, label) {
  data.frame(
    cohort = label,
    n = nrow(dat),
    events = sum(dat$status),
    median_time = round(median(dat$time), 2),
    mean_X1 = round(mean(dat$X_1), 2),
    mean_X2 = round(mean(dat$X_2), 2),
    mean_X3 = round(mean(dat$X_3), 2)
  )
}
sample_tbl <- rbind(
  cohort_line(historical, "Historical controls"),
  cohort_line(trial[trial$arm == 0L, ], "Concurrent controls"),
  cohort_line(trial[trial$arm == 1L, ], "Treated")
)

# Pooled unadjusted
keep <- c("id", "X_1", "X_2", "X_3", "X_4", "X_5", "X_6", "time", "status", "trt", "s.ob")
historical$trt <- 0
pooled <- rbind(historical[, keep], trial[, keep])
pooled$s.ob <- survival::Surv(pooled$time, pooled$status)
pooled_fit <- survival::coxph(s.ob ~ trt, data = pooled)
pooled_co <- summary(pooled_fit)$coefficients
pooled_est <- unname(pooled_co["trt", "coef"])
pooled_se <- unname(pooled_co["trt", "se(coef)"])
pooled_z <- stats::qnorm(0.975)

# Population synthetic control
cov_keep <- c("X_1", "X_2", "X_3", "X_4", "X_5", "X_6", "time", "status")
bal <- entropy_balance_weights(historical[, cov_keep], trial[, cov_keep])
historical$w <- bal$w
trial$w <- 1
sc_keep <- c(keep, "w")
sc <- rbind(historical[, sc_keep], trial[, sc_keep])
sc$s.ob <- survival::Surv(sc$time, sc$status)
sc_fit <- survival::coxph(s.ob ~ trt, data = sc, weights = w)
sc_co <- summary(sc_fit)$coefficients
sc_est <- unname(sc_co["trt", "coef"])
sc_se <- unname(sc_co["trt", "se(coef)"])

# Bayesian case-weighted
set.seed(seed_h)
a0 <- bal$w
log_post_cw <- function(theta) {
  lambda0 <- theta[1]
  beta <- theta[2]
  lp <- stats::dnorm(lambda0, 0, 10, log = TRUE) + stats::dnorm(beta, 0, 10, log = TRUE)
  rate_c <- exp(-lambda0)
  rate_t <- exp(-lambda0 + beta)
  conc <- trial$arm == 0L
  trt <- trial$arm == 1L
  lp <- lp + sum(trial$status[conc] * log(rate_c) - rate_c * trial$time[conc])
  lp <- lp + sum(trial$status[trt] * log(rate_t) - rate_t * trial$time[trt])
  lp <- lp + sum(a0 * (historical$status * log(rate_c) - rate_c * historical$time))
  lp
}
theta <- c(0, 0)
cur <- log_post_cw(theta)
step <- c(0.12, 0.25)
draws_cw <- matrix(NA_real_, 8000L, 2L)
acc <- 0L
for (i in seq_len(8000L)) {
  prop <- theta + stats::rnorm(2, 0, step)
  pr <- log_post_cw(prop)
  if (is.finite(pr) && log(stats::runif(1)) < (pr - cur)) {
    theta <- prop
    cur <- pr
    acc <- acc + 1L
  }
  draws_cw[i, ] <- theta
  if (i %% 500L == 0L && i < 2000L) {
    rate <- acc / 500
    if (rate < 0.15) step <- step * 0.7
    if (rate > 0.45) step <- step * 1.3
    acc <- 0L
  }
}
beta_cw <- draws_cw[2001:8000, 2]
cw_q <- stats::quantile(beta_cw, c(0.025, 0.975))

# Commensurate prior, covariate-adjusted
sr <- survival::survreg(
  survival::Surv(time, status) ~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6,
  data = historical,
  dist = "exponential"
)
mu_hist <- unname(stats::coef(sr)["(Intercept)"])
X_trial <- stats::model.matrix(~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6, data = trial)[, -1, drop = FALSE]
X_hist <- stats::model.matrix(~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6, data = historical)[, -1, drop = FALSE]
set.seed(seed_h + 3L)
log_post_comm <- function(theta) {
  lambda0 <- theta[1]
  gamma <- theta[2:10]
  beta <- theta[11]
  tau <- exp(theta[12])
  lp <- stats::dnorm(lambda0, mu_hist, 1 / sqrt(tau), log = TRUE)
  lp <- lp + sum(stats::dnorm(gamma, 0, 10, log = TRUE))
  lp <- lp + stats::dnorm(beta, 0, 10, log = TRUE)
  lp <- lp + stats::dgamma(tau, 1, 1, log = TRUE) + theta[12]
  eta <- -lambda0 + as.numeric(X_trial %*% gamma) + beta * trial$arm
  lp <- lp + sum(trial$status * eta - exp(eta) * trial$time)
  lp
}
theta <- c(mu_hist, -unname(stats::coef(sr)[-1]), log(0.7), 0)
cur <- log_post_comm(theta)
step <- rep(0.08, 12L)
draws_comm <- matrix(NA_real_, 6000L, 12L)
acc <- 0L
for (i in seq_len(6000L)) {
  prop <- theta + stats::rnorm(12L, 0, step)
  pr <- log_post_comm(prop)
  if (is.finite(pr) && log(stats::runif(1)) < (pr - cur)) {
    theta <- prop
    cur <- pr
    acc <- acc + 1L
  }
  draws_comm[i, ] <- theta
  if (i %% 400L == 0L && i < 1500L) {
    rate <- acc / 400
    if (rate < 0.15) step <- step * 0.75
    if (rate > 0.35) step <- step * 1.15
    acc <- 0L
  }
}
keep_c <- draws_comm[1501:6000, , drop = FALSE]
comm_beta <- keep_c[, 11]
comm_tau <- exp(keep_c[, 12])
comm_l0 <- keep_c[, 1]

# LEAP, two classes, covariate-adjusted shared hazard
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
beta_leap <- rep(NA_real_, n_leap)
g1_leap <- rep(NA_real_, n_leap)
n1_leap <- rep(NA_integer_, n_leap)
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
    # class-2 likelihood plus prior; theta part unchanged
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
  beta_leap[i] <- theta_l[11]
  g1_leap[i] <- gamma1
  n1_leap[i] <- n1
  if (i %% 200L == 0L && i < burn_l) {
    rate <- acc_l / 200
    if (rate < 0.15) step_l <- step_l * 0.75
    if (rate > 0.35) step_l <- step_l * 1.1
    acc_l <- 0L
  }
}
beta_l <- beta_leap[(burn_l + 1L):n_leap]
g1_l <- g1_leap[(burn_l + 1L):n_leap]
p_bar <- colMeans(p_exch[(burn_l + 1L):n_leap, , drop = FALSE])
ord <- order(p_bar)
leap_pts <- data.frame(
  rank = c(1:3, 498:500),
  p_exchangeable = round(p_bar[c(ord[1:3], ord[498:500])], 3),
  time = round(historical$time[c(ord[1:3], ord[498:500])], 2),
  event = historical$status[c(ord[1:3], ord[498:500])],
  X_1 = round(historical$X_1[c(ord[1:3], ord[498:500])], 2)
)

# PSC: calibrate, then adjusted joint weights, then combine the two contrasts
set.seed(seed_h + 1L)
train_id <- sample.int(500L, floor(0.75 * 500))
train <- historical[train_id, ]
valid <- historical[-train_id, ]
valid$cen <- valid$status
train_fit <- flexsurv::flexsurvspline(s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6, data = train, k = 1)
cal_fit <- psc::pscfit(as_pscCFM_dev(train_fit), valid, nsim = 2000L, nchain = 1L)
cal_beta <- as.numeric(posterior::as_draws_df(cal_fit$draws)$beta_1)
cal_beta <- cal_beta[is.finite(cal_beta)]
cal_hr <- stats::median(exp(cal_beta))
psc_ok <- is.finite(cal_hr) && cal_hr >= 0.8 && cal_hr <= 1.2
psc_out <- list(cal_hr = cal_hr, accepted = psc_ok, n_train = nrow(train), n_valid = nrow(valid),
                events_train = sum(train$status), events_valid = sum(valid$status))
if (psc_ok) {
  cfm <- flexsurv::flexsurvspline(s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6, data = historical, k = 1)
  cfmw <- as_pscCFM_dev(cfm)
  pw <- pscSupportWeights(cfmw, trial, method = "joint")
  pw$weights <- pmin(1, pw$weights / 0.5)
  trt_fac <- factor(ifelse(trial$arm == 0L, "control", "treated"), levels = c("control", "treated"))
  set.seed(seed_h + 2L)
  psc_fit <- psc::pscfit(cfmw, trial, nsim = 2000L, nchain = 1L, proximity_weights = pw, trt = trt_fac)
  drs <- posterior::as_draws_df(psc_fit$draws)
  b1 <- as.numeric(drs$beta_1)
  b2 <- as.numeric(drs$beta_2)
  okb <- is.finite(b1) & is.finite(b2)
  b1 <- b1[okb]
  b2 <- b2[okb]
  direct <- b2
  indirect <- b2 - b1
  sig_d <- stats::var(direct)
  sig_i <- stats::var(indirect)
  rho <- stats::cor(direct, indirect)
  w_d <- sig_i / (sig_d + sig_i)
  mu_comb <- w_d * mean(direct) + (1 - w_d) * mean(indirect)
  sig_comb <- w_d^2 * sig_d + (1 - w_d)^2 * sig_i + 2 * w_d * (1 - w_d) * rho * sqrt(sig_d) * sqrt(sig_i)
  psc_out$direct <- mean(direct)
  psc_out$indirect <- mean(indirect)
  psc_out$weight_direct <- w_d
  psc_out$combined <- mu_comb
  psc_out$combined_se <- sqrt(sig_comb)
  psc_out$sum_w <- sum(pw$weights)
  psc_out$mean_w_control <- mean(pw$weights[trial$arm == 0L])
  psc_out$mean_w_treated <- mean(pw$weights[trial$arm == 1L])
  psc_out$beta_names <- grep("^beta_", names(drs), value = TRUE)
}

z <- stats::qnorm(0.975)
if (!psc_ok) {
  psc_out$direct <- NA_real_
  psc_out$indirect <- NA_real_
  psc_out$weight_direct <- NA_real_
  psc_out$combined <- NA_real_
  psc_out$combined_se <- NA_real_
}
data.frame(
  est_pooled = pooled_est,
  lo_pooled = pooled_est - z * pooled_se,
  hi_pooled = pooled_est + z * pooled_se,
  est_sc = sc_est,
  lo_sc = sc_est - z * sc_se,
  hi_sc = sc_est + z * sc_se,
  ess_sc = bal$ess,
  est_cw = stats::median(beta_cw),
  lo_cw = unname(cw_q[1]),
  hi_cw = unname(cw_q[2]),
  est_psc = psc_out$combined,
  lo_psc = psc_out$combined - z * psc_out$combined_se,
  hi_psc = psc_out$combined + z * psc_out$combined_se,
  psc_accepted = psc_ok,
  cal_hr = cal_hr,
  psc_direct = psc_out$direct,
  psc_indirect = psc_out$indirect,
  psc_w_direct = psc_out$weight_direct,
  est_comm = stats::median(comm_beta),
  lo_comm = unname(stats::quantile(comm_beta, 0.025)),
  hi_comm = unname(stats::quantile(comm_beta, 0.975)),
  tau = stats::median(comm_tau),
  est_leap = stats::median(beta_l),
  lo_leap = unname(stats::quantile(beta_l, 0.025)),
  hi_leap = unname(stats::quantile(beta_l, 0.975)),
  gamma1 = stats::median(g1_l),
  n_exchangeable = stats::median(n1_leap[(burn_l + 1L):n_leap])
)
}

truth <- log(0.7)
rep_from <- 11L
rep_to <- 20L
rows <- vector("list", rep_to - rep_from + 1L)
for (rep_i in rep_from:rep_to) {
  seed_h <- 11L + (rep_i - 1L) * 1000L
  rows[[rep_i - rep_from + 1L]] <- run_rep(seed_h)
  rows[[rep_i - rep_from + 1L]]$rep <- rep_i
  cat("Finished", rep_i, "\n")
  flush.console()
}
new_rows <- do.call(rbind, rows)
prev_path <- file.path(sim_dir, "data/hybrid_10.csv")
tab <- rbind(utils::read.csv(prev_path), new_rows)
out_path <- file.path(sim_dir, "data/hybrid_20.csv")
utils::write.csv(tab, out_path, row.names = FALSE)
cover <- function(est, lo, hi) {
  use <- is.finite(est)
  c(
    n = sum(use),
    mean = mean(est[use]),
    median = stats::median(est[use]),
    bias = mean(est[use] - truth),
    cover = sum(lo[use] <= truth & hi[use] >= truth)
  )
}
summary_tbl <- rbind(
  pooled = cover(tab$est_pooled, tab$lo_pooled, tab$hi_pooled),
  population_sc = cover(tab$est_sc, tab$lo_sc, tab$hi_sc),
  case_weighted = cover(tab$est_cw, tab$lo_cw, tab$hi_cw),
  psc = cover(tab$est_psc, tab$lo_psc, tab$hi_psc),
  commensurate = cover(tab$est_comm, tab$lo_comm, tab$hi_comm),
  leap = cover(tab$est_leap, tab$lo_leap, tab$hi_leap)
)
cat("Truth", truth, "PSC accepted", sum(tab$psc_accepted), "of", nrow(tab), "\n")
print(round(summary_tbl, 3))
cat("Saved", out_path, "\n")
