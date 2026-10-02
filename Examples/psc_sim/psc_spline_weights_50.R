# Same-margin samples on the spline master.
# Historical controls and the treated trial are both drawn from dataset 1.
# PSC is fit only when the validation hazard ratio is inside 0.8 to 1.2.
# Each accepted sample is fit with three proximity weights, from the same trial seed:
# marginal geometric mean, joint Mahalanobis, and joint inflated by min(1, w/0.5).

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

sim_dir <- "~/Documents/GitHub/synthetic-Controlled-Trials/Examples/psc_sim"
sim_dir <- normalizePath(sim_dir, mustWork = TRUE)
source(file.path(sim_dir, "approaches/psc.R"))
load(file.path(sim_dir, "data/master.RData"))

n_hist <- 500L
n_trial <- 75L
rep_from <- 151L
rep_to <- 350L
n_cores <- 4L
nsim <- 2000L
true_log_hr <- log(0.7)

prepare_arm <- function(dat, outcome = c("control", "treated")) {
  outcome <- match.arg(outcome)
  if (outcome == "control") {
    dat$time <- dat$t0
    dat$status <- dat$event0
  } else {
    dat$time <- dat$t1
    dat$status <- dat$event1
  }
  dat
}

sample_cohort <- function(dat, n, seed) {
  set.seed(seed)
  dat[sample.int(nrow(dat), n), , drop = FALSE]
}

load_psc_for_prototype()
assignInNamespace(
  "txtProgressBar",
  function(...) structure(list(), class = "txtProgressBar"),
  "utils"
)
assignInNamespace(
  "setTxtProgressBar",
  function(pb, ...) invisible(pb),
  "utils"
)

beta_summary <- function(fit, weights) {
  beta <- as.numeric(posterior::as_draws_df(fit$draws)$beta_1)
  beta <- beta[is.finite(beta)]
  q <- stats::quantile(beta, c(0.025, 0.975), names = FALSE)
  w <- weights[is.finite(weights)]
  c(est = stats::median(beta), lo = unname(q[1]), hi = unname(q[2]), sum_w = sum(w))
}

run_one <- function(rep_i) {
  seed_h <- 11L + (rep_i - 1L) * 1000L
  out <- data.frame(
    rep = rep_i,
    accepted = FALSE,
    calibration_hr = NA_real_,
    est_g = NA_real_, lo_g = NA_real_, hi_g = NA_real_, sum_g = NA_real_,
    est_j = NA_real_, lo_j = NA_real_, hi_j = NA_real_, sum_j = NA_real_,
    est_i = NA_real_, lo_i = NA_real_, hi_i = NA_real_, sum_i = NA_real_
  )
  utils::capture.output({
    historical <- prepare_arm(
      sample_cohort(master$dataset_1, n_hist, seed_h),
      "control"
    )
    trial <- prepare_arm(
      sample_cohort(master$dataset_1, n_trial, seed_h + 100L),
      "treated"
    )
    historical$s.ob <- survival::Surv(historical$time, historical$status)
    historical$X_4 <- factor(historical$X_4, levels = 1:3)
    historical$X_5 <- factor(historical$X_5, levels = 1:3)
    historical$X_6 <- factor(historical$X_6, levels = 1:3)
    trial$X_4 <- factor(trial$X_4, levels = levels(historical$X_4))
    trial$X_5 <- factor(trial$X_5, levels = levels(historical$X_5))
    trial$X_6 <- factor(trial$X_6, levels = levels(historical$X_6))
    trial$cen <- trial$status
    trial$s.ob <- survival::Surv(trial$time, trial$cen)

    set.seed(seed_h + 1L)
    train_id <- sample.int(n_hist, floor(0.75 * n_hist))
    train <- historical[train_id, , drop = FALSE]
    valid <- historical[-train_id, , drop = FALSE]
    train_fit <- flexsurv::flexsurvspline(
      s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6,
      data = train,
      k = 1
    )
    valid$cen <- valid$status
    valid$s.ob <- survival::Surv(valid$time, valid$cen)
    cal_fit <- psc::pscfit(
      as_pscCFM_dev(train_fit),
      valid,
      nsim = nsim,
      nchain = 1L
    )
    cal_beta <- as.numeric(posterior::as_draws_df(cal_fit$draws)$beta_1)
    cal_beta <- cal_beta[is.finite(cal_beta)]
    cal_hr <- if (length(cal_beta) == 0L) NA_real_ else stats::median(exp(cal_beta))
    out$calibration_hr <- cal_hr
    accepted <- is.finite(cal_hr) && cal_hr >= 0.8 && cal_hr <= 1.2
    out$accepted <- accepted
    if (accepted) {
      cfm <- flexsurv::flexsurvspline(
        s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6,
        data = historical,
        k = 1
      )
      cfmw <- as_pscCFM_dev(cfm)
      w_all <- pscSupportWeights(cfmw, trial, method = "all")
      pw_g <- w_all$marginal_geomean
      pw_j <- w_all$joint_cdf
      pw_i <- pw_j
      pw_i$weights <- pmin(1, pw_j$weights / 0.5)
      pw_i$method <- "joint_inflated"

      set.seed(seed_h + 2L)
      fit_g <- psc::pscfit(cfmw, trial, nsim = nsim, nchain = 1L, proximity_weights = pw_g)
      set.seed(seed_h + 2L)
      fit_j <- psc::pscfit(cfmw, trial, nsim = nsim, nchain = 1L, proximity_weights = pw_j)
      set.seed(seed_h + 2L)
      fit_i <- psc::pscfit(cfmw, trial, nsim = nsim, nchain = 1L, proximity_weights = pw_i)

      sg <- beta_summary(fit_g, pw_g$weights)
      sj <- beta_summary(fit_j, pw_j$weights)
      si <- beta_summary(fit_i, pw_i$weights)
      out$est_g <- sg["est"]; out$lo_g <- sg["lo"]; out$hi_g <- sg["hi"]; out$sum_g <- sg["sum_w"]
      out$est_j <- sj["est"]; out$lo_j <- sj["lo"]; out$hi_j <- sj["hi"]; out$sum_j <- sj["sum_w"]
      out$est_i <- si["est"]; out$lo_i <- si["lo"]; out$hi_i <- si["hi"]; out$sum_i <- si["sum_w"]
    }
  })
  out
}

reps <- seq.int(rep_from, rep_to)
out_path <- file.path(sim_dir, "data/psc_spline_three_weights_151_350.csv")
results <- vector("list", length(reps))
batch <- n_cores * 2L
for (start in seq_along(reps)) {
  if ((start - 1L) %% batch != 0L) next
  part <- start:min(start + batch - 1L, length(reps))
  t0 <- proc.time()[["elapsed"]]
  results[part] <- parallel::mclapply(reps[part], run_one, mc.cores = n_cores)
  tab <- do.call(rbind, results[seq_len(max(part))])
  utils::write.csv(tab, out_path, row.names = FALSE)
  acc <- sum(tab$accepted, na.rm = TRUE)
  cat(
    "Finished rep", max(reps[part]),
    " done", max(part), "of", length(reps),
    " accepted", acc,
    " elapsed_batch_sec", round(proc.time()[["elapsed"]] - t0, 1),
    "\n"
  )
  flush.console()
}

tab <- do.call(rbind, results)
utils::write.csv(tab, out_path, row.names = FALSE)
first <- utils::read.csv(file.path(sim_dir, "data/psc_spline_three_weights_150.csv"))
tab <- rbind(first, tab)
out_path <- file.path(sim_dir, "data/psc_spline_three_weights_350.csv")
utils::write.csv(tab, out_path, row.names = FALSE)
ok <- tab[tab$accepted, ]
report <- function(est, lo, hi, sum_w) {
  data.frame(
    mean = mean(est),
    median = stats::median(est),
    sd = stats::sd(est),
    mean_bias = mean(est - true_log_hr),
    coverage = sum(lo <= true_log_hr & hi >= true_log_hr),
    mean_sum_w = mean(sum_w),
    min_sum_w = min(sum_w),
    max_sum_w = max(sum_w)
  )
}
cat("Accepted", nrow(ok), "of", nrow(tab), " truth", true_log_hr, "\n")
print(rbind(
  geometric_mean = report(ok$est_g, ok$lo_g, ok$hi_g, ok$sum_g),
  joint = report(ok$est_j, ok$lo_j, ok$hi_j, ok$sum_j),
  joint_inflated = report(ok$est_i, ok$lo_i, ok$hi_i, ok$sum_i)
))
cat("Saved", out_path, "\n")
