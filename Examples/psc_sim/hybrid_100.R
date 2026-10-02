# Hybrid simulation across 4 cores.
# PSC uses the corrected treatment effect beta_2 - lambda * beta_1.
# Commensurate prior and LEAP use a one-knot flexsurvspline baseline.
# 100 replicates of every method on all five trial populations.
# 500 historical controls from dataset 1.
# 75 trial patients, randomized 38:37.
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
try(assignInNamespace("txtProgressBar", function(...) structure(list(), class = "txtProgressBar"), "utils"), silent = TRUE)
try(assignInNamespace("setTxtProgressBar", function(pb, ...) invisible(pb), "utils"), silent = TRUE)

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

run_rep <- function(seed_h, trial_name) {
set.seed(seed_h)
historical <- arm_of(master$dataset_1[sample.int(nrow(master$dataset_1), 500L), ], "control")
set.seed(seed_h + 100L)
trial_source <- master[[trial_name]]
trial_rows <- trial_source[sample.int(nrow(trial_source), 75L), ]
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

# One-knot Royston-Parmar baseline, the flexsurvspline model on the hazard scale.
# Knots are the ones chosen by flexsurv from the historical controls.
# gamma0, gamma1 and gamma2 are the baseline. Covariates add to gamma0.
hist_fit <- flexsurv::flexsurvspline(
  s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6,
  data = historical,
  k = 1
)
spline_knots <- hist_fit$knots
coef_hist <- unname(stats::coef(hist_fit))
base_hist <- coef_hist[1:3]
cov_f <- ~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6
X_hist <- stats::model.matrix(cov_f, data = historical)[, -1L, drop = FALSE]
X_trial <- stats::model.matrix(cov_f, data = trial)[, -1L, drop = FALSE]
spline_ll_i <- function(time, status, gamma0, shape, lp) {
  gamma <- cbind(gamma0 + lp, shape[1], shape[2])
  ld <- flexsurv::dsurvspline(time, gamma = gamma, knots = spline_knots, scale = "hazard", log = TRUE)
  ls <- flexsurv::psurvspline(
    time, gamma = gamma, knots = spline_knots, scale = "hazard",
    lower.tail = FALSE, log.p = TRUE
  )
  ll <- ifelse(status == 1, ld, ls)
  ll[!is.finite(ll)] <- -1e8
  ll
}

# Commensurate prior on the three baseline spline coefficients.
# Prior mean is the historical flexsurvspline fit. Prior SD is 1/sqrt(tau).
# The trial likelihood uses that spline. Historical patients do not re-enter it.
set.seed(seed_h + 3L)
log_post_comm <- function(theta) {
  base <- theta[1:3]
  gx <- theta[4:12]
  beta <- theta[13]
  tau <- exp(theta[14])
  lp <- sum(stats::dnorm(base, base_hist, 1 / sqrt(tau), log = TRUE))
  lp <- lp + sum(stats::dnorm(gx, 0, 10, log = TRUE))
  lp <- lp + stats::dnorm(beta, 0, 10, log = TRUE)
  lp <- lp + stats::dgamma(tau, 1, 1, log = TRUE) + theta[14]
  eta_lp <- as.numeric(X_trial %*% gx) + beta * trial$arm
  lp + sum(spline_ll_i(trial$time, trial$status, base[1], base[2:3], eta_lp))
}
theta <- c(coef_hist, log(0.7), 0)
cur <- log_post_comm(theta)
step <- rep(0.05, 14L)
draws_comm <- matrix(NA_real_, 6000L, 14L)
acc <- 0L
for (i in seq_len(6000L)) {
  prop <- theta + stats::rnorm(14L, 0, step)
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
comm_beta <- keep_c[, 13]
comm_tau <- exp(keep_c[, 14])

# LEAP. The exchangeable class shares the trial spline, including covariates.
# The other class has its own spline baseline and no covariates.
set.seed(seed_h + 4L)
n_h <- nrow(historical)
gamma1 <- 0.8
base2 <- base_hist
theta_l <- c(coef_hist, log(0.7))
loglik_class <- function(theta, base2) {
  ll1 <- spline_ll_i(
    historical$time, historical$status, theta[1], theta[2:3],
    as.numeric(X_hist %*% theta[4:12])
  )
  ll2 <- spline_ll_i(historical$time, historical$status, base2[1], base2[2:3], 0)
  cbind(ll1, ll2)
}
log_post_leap <- function(theta, class1, base2) {
  lp <- sum(stats::dnorm(theta[1:13], 0, 10, log = TRUE))
  lp <- lp + sum(stats::dnorm(base2, 0, 10, log = TRUE))
  eta_lp <- as.numeric(X_trial %*% theta[4:12]) + theta[13] * trial$arm
  lp <- lp + sum(spline_ll_i(trial$time, trial$status, theta[1], theta[2:3], eta_lp))
  if (any(class1)) {
    lp <- lp + sum(spline_ll_i(
      historical$time[class1], historical$status[class1], theta[1], theta[2:3],
      as.numeric(X_hist[class1, , drop = FALSE] %*% theta[4:12])
    ))
  }
  lp
}
class1 <- rep(TRUE, n_h)
cur <- log_post_leap(theta_l, class1, base2)
step_l <- rep(0.05, 13L)
step2 <- 0.05
n_leap <- 4000L
burn_l <- 1000L
p_exch <- matrix(NA_real_, n_leap, n_h)
beta_leap <- rep(NA_real_, n_leap)
g1_leap <- rep(NA_real_, n_leap)
n1_leap <- rep(NA_integer_, n_leap)
acc_l <- 0L
for (i in seq_len(n_leap)) {
  ll <- loglik_class(theta_l, base2)
  logp1 <- log(gamma1) + ll[, 1]
  logp2 <- log(1 - gamma1) + ll[, 2]
  m <- pmax(logp1, logp2)
  p1 <- exp(logp1 - m) / (exp(logp1 - m) + exp(logp2 - m))
  class1 <- stats::runif(n_h) < p1
  p_exch[i, ] <- p1
  n1 <- sum(class1)
  gamma1 <- stats::rbeta(1, 1 + n1, 1 + (n_h - n1))
  cur <- log_post_leap(theta_l, class1, base2)
  prop <- theta_l + stats::rnorm(13L, 0, step_l)
  pr <- log_post_leap(prop, class1, base2)
  if (is.finite(pr) && log(stats::runif(1)) < (pr - cur)) {
    theta_l <- prop
    cur <- pr
    acc_l <- acc_l + 1L
  }
  if (any(!class1)) {
    prop2 <- base2 + stats::rnorm(3L, 0, step2)
    use <- !class1
    ll_old <- sum(spline_ll_i(historical$time[use], historical$status[use], base2[1], base2[2:3], 0)) +
      sum(stats::dnorm(base2, 0, 10, log = TRUE))
    ll_new <- sum(spline_ll_i(historical$time[use], historical$status[use], prop2[1], prop2[2:3], 0)) +
      sum(stats::dnorm(prop2, 0, 10, log = TRUE))
    if (is.finite(ll_new) && log(stats::runif(1)) < (ll_new - ll_old)) {
      base2 <- prop2
      cur <- log_post_leap(theta_l, class1, base2)
    }
  }
  beta_leap[i] <- theta_l[13]
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
  # beta_2 is treated versus the historical model.
  # beta_2 - beta_1 is treated versus the concurrent controls.
  # The reported treatment effect subtracts only the part of beta_1 that
  # exceeds its posterior noise: beta_2 - lambda * beta_1.
  historical <- b2
  within_trial <- b2 - b1
  m1 <- mean(b1)
  v1 <- stats::var(b1)
  lambda <- 0
  if (is.finite(m1) && is.finite(v1) && m1^2 > v1) lambda <- 1 - v1 / m1^2
  corrected <- b2 - lambda * b1
  psc_out$direct <- mean(within_trial)
  psc_out$indirect <- mean(historical)
  psc_out$beta1 <- m1
  psc_out$lambda <- lambda
  psc_out$combined <- mean(corrected)
  psc_out$combined_se <- stats::sd(corrected)
  psc_out$sum_w <- sum(pw$weights)
  psc_out$mean_w_control <- mean(pw$weights[trial$arm == 0L])
  psc_out$mean_w_treated <- mean(pw$weights[trial$arm == 1L])
  psc_out$beta_names <- grep("^beta_", names(drs), value = TRUE)
}

z <- stats::qnorm(0.975)
if (!psc_ok) {
  psc_out$direct <- NA_real_
  psc_out$indirect <- NA_real_
  psc_out$combined <- NA_real_
  psc_out$combined_se <- NA_real_
  psc_out$beta1 <- NA_real_
  psc_out$lambda <- NA_real_
  psc_out$sum_w <- NA_real_
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
  psc_beta1 = psc_out$beta1,
  psc_lambda = psc_out$lambda,
  sum_w = psc_out$sum_w,
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

if (!identical(Sys.getenv("HYBRID_SOURCE_ONLY"), "1")) {
truth <- log(0.7)
n_reps <- 100L
n_cores <- 4L
batch_reps <- 10L
overlaps <- c(
  dataset_1 = "Same margins",
  dataset_2 = "Moderate overlap",
  dataset_3 = "Poor overlap",
  dataset_2_shift = "Moderate overlap, outcome shift",
  dataset_3_shift = "Poor overlap, outcome shift"
)
seeds <- c(dataset_1 = 11L, dataset_2 = 12L, dataset_3 = 13L, dataset_2_shift = 14L, dataset_3_shift = 15L)
out_path <- file.path(sim_dir, "data/hybrid_100.csv")
done_df <- NULL

run_one <- function(i) {
  rep_i <- tasks$rep_i[i]
  nm <- tasks$nm[i]
  seed_h <- unname(seeds[[nm]]) + (rep_i - 1L) * 1000L
  one <- NULL
  utils::capture.output({
    one <- run_rep(seed_h, nm)
  })
  one$rep <- rep_i
  one$overlap <- unname(overlaps[[nm]])
  one
}

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
n_ready <- function(sub) {
  min(
    sum(is.finite(sub$est_pooled)),
    sum(is.finite(sub$est_sc)),
    sum(is.finite(sub$est_cw)),
    sum(is.finite(sub$est_psc)),
    sum(is.finite(sub$est_comm)),
    sum(is.finite(sub$est_leap))
  )
}
report <- function(tab) {
  for (ov in overlaps) {
    sub <- tab[tab$overlap == ov, ]
    summary_tbl <- rbind(
      pooled = cover(sub$est_pooled, sub$lo_pooled, sub$hi_pooled),
      population_sc = cover(sub$est_sc, sub$lo_sc, sub$hi_sc),
      case_weighted = cover(sub$est_cw, sub$lo_cw, sub$hi_cw),
      psc = cover(sub$est_psc, sub$lo_psc, sub$hi_psc),
      commensurate = cover(sub$est_comm, sub$lo_comm, sub$hi_comm),
      leap = cover(sub$est_leap, sub$lo_leap, sub$hi_leap)
    )
    cat(ov, "replicates", nrow(sub), "PSC accepted", sum(sub$psc_accepted),
        "mean lambda", round(mean(sub$psc_lambda, na.rm = TRUE), 3), "\n")
    print(round(summary_tbl, 3))
  }
}

next_rep <- 1L
while (next_rep <= n_reps) {
  rep_to <- min(next_rep + batch_reps - 1L, n_reps)
  tasks <- expand.grid(
    rep_i = next_rep:rep_to,
    nm = names(overlaps),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  t0 <- proc.time()[["elapsed"]]
  new_rows <- parallel::mclapply(seq_len(nrow(tasks)), run_one, mc.cores = n_cores)
  new_df <- do.call(rbind, new_rows)
  done_df <- if (is.null(done_df)) new_df else rbind(done_df, new_df[, names(done_df), drop = FALSE])
  utils::write.csv(done_df, out_path, row.names = FALSE)
  cat("UPDATE through replicate", rep_to, "sec", round(proc.time()[["elapsed"]] - t0, 1), "\n")
  report(done_df)
  flush.console()
  next_rep <- rep_to + 1L
}
cat("Truth", truth, "\n")
cat("Saved", out_path, "\n")
}
