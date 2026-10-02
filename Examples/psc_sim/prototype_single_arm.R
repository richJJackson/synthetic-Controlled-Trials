# Single-arm simulation: 500 replicates at each overlap, across 4 cores.
# PSC uses adjusted joint proximity weights, min(1, w/0.5).
# PSC is fit only when the validation hazard ratio is inside 0.8 to 1.2.
# Historical controls are sampled from dataset 1.
# The trial sample is taken from dataset 1, 2, or 3 and all receive treatment.

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

sim_dir <- "~/Documents/GitHub/synthetic-Controlled-Trials/Examples/psc_sim"
sim_dir <- normalizePath(sim_dir, mustWork = TRUE)
source(file.path(sim_dir, "approaches/entropy_balance.R"))
source(file.path(sim_dir, "approaches/pooled_unadjusted.R"))
source(file.path(sim_dir, "approaches/population_sc.R"))
source(file.path(sim_dir, "approaches/bayes_case_weighted.R"))
source(file.path(sim_dir, "approaches/psc.R"))

load(file.path(sim_dir, "data/master.RData"))

n_hist <- 500L
n_trial <- 75L
rep_from <- 1L
rep_to <- 500L
n_cores <- 4L
true_log_hr <- log(0.7)
overlaps <- c(dataset_1 = "Same margins", dataset_2 = "Moderate overlap", dataset_3 = "Poor overlap")
seeds <- c(dataset_1 = 11L, dataset_2 = 12L, dataset_3 = 13L)

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

tasks <- expand.grid(
  rep_i = seq.int(rep_from, rep_to),
  nm = names(overlaps),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

run_one <- function(i) {
  rep_i <- tasks$rep_i[i]
  nm <- tasks$nm[i]
  seed_h <- unname(seeds[[nm]]) + (rep_i - 1L) * 1000L
  one <- NULL
  utils::capture.output({
    historical <- prepare_arm(
      sample_cohort(master$dataset_1, n_hist, seed_h),
      "control"
    )
    trial <- prepare_arm(
      sample_cohort(master[[nm]], n_trial, seed_h + 100L),
      "treated"
    )
    one <- rbind(
      fit_pooled_unadjusted(historical, trial),
      fit_population_sc(historical, trial),
      fit_bayes_case_weighted(historical, trial, seed = seed_h),
      fit_psc(historical, trial, seed = seed_h + 1L)
    )
  })
  one$rep <- rep_i
  one$overlap <- overlaps[[nm]]
  one$trial_dataset <- nm
  one$true_log_hr <- true_log_hr
  one$bias <- one$est - true_log_hr
  one$covered <- one$lo <= true_log_hr & one$hi >= true_log_hr
  one
}

batch <- 50L
idx <- seq_len(nrow(tasks))
results <- vector("list", length(idx))
reps_out <- file.path(sim_dir, "data/single_arm_adjusted_joint_reps.csv")
for (start in seq(1L, length(idx), by = batch)) {
  part <- idx[start:min(start + batch - 1L, length(idx))]
  t0 <- proc.time()[["elapsed"]]
  results[part] <- parallel::mclapply(part, run_one, mc.cores = n_cores)
  done <- do.call(rbind, results[seq_len(max(part))])
  utils::write.csv(done, reps_out, row.names = FALSE)
  cat("UPDATE", max(part), "of", length(idx), "sec", round(proc.time()[["elapsed"]] - t0, 1), "\n")
  show <- done
  show$covered_n <- as.integer(show$lo <= show$true_log_hr & show$hi >= show$true_log_hr)
  for (ov in unique(show$overlap)) {
    for (me in unique(show$method)) {
      sub <- show[show$overlap == ov & show$method == me, ]
      use <- if (me == "PSC proximity-weighted") sub$accepted %in% TRUE else rep(TRUE, nrow(sub))
      est <- sub$est[use]
      cat(
        ov, "|", me,
        "| n", sum(use),
        "| mean", round(mean(est), 3),
        "| bias", round(mean(est - true_log_hr), 3),
        "| cover", sum(sub$covered_n[use], na.rm = TRUE), "/", sum(use),
        "| sumw", round(mean(sub$sum_w[use], na.rm = TRUE), 1),
        "\n"
      )
    }
  }
  flush.console()
}

prototype <- do.call(rbind, results)
rownames(prototype) <- NULL
reps_out <- file.path(sim_dir, "data/prototype_single_arm_reps.csv")
prototype$overlap <- factor(prototype$overlap, levels = overlaps)
prototype$method <- factor(prototype$method, levels = c(
  "Pooled unadjusted",
  "Population synthetic control",
  "Bayesian case-weighted",
  "PSC proximity-weighted"
))

mean_ok <- function(x) mean(x, na.rm = TRUE)
means <- aggregate(
  cbind(est, se, bias, sum_w, ess, covered) ~ overlap + method,
  data = prototype,
  FUN = mean_ok
)
sds <- aggregate(est ~ overlap + method, data = prototype, FUN = function(x) stats::sd(x, na.rm = TRUE))
names(sds)[names(sds) == "est"] <- "mc_sd"
med_est <- aggregate(est ~ overlap + method, data = prototype, FUN = function(x) stats::median(x, na.rm = TRUE))
names(med_est)[names(med_est) == "est"] <- "median_est"
med_ess <- aggregate(ess ~ overlap + method, data = prototype, FUN = function(x) stats::median(x, na.rm = TRUE))
names(med_ess)[names(med_ess) == "ess"] <- "median_ess"
summary_tbl <- merge(means, sds, by = c("overlap", "method"))
summary_tbl <- merge(summary_tbl, med_est, by = c("overlap", "method"))
summary_tbl <- merge(summary_tbl, med_ess, by = c("overlap", "method"))
summary_tbl <- summary_tbl[order(summary_tbl$overlap, summary_tbl$method), ]
names(summary_tbl)[names(summary_tbl) == "est"] <- "mean_est"
names(summary_tbl)[names(summary_tbl) == "se"] <- "mean_se"
names(summary_tbl)[names(summary_tbl) == "bias"] <- "mean_bias"
names(summary_tbl)[names(summary_tbl) == "sum_w"] <- "mean_sum_w"
names(summary_tbl)[names(summary_tbl) == "ess"] <- "mean_ess"
names(summary_tbl)[names(summary_tbl) == "covered"] <- "coverage"
summary_tbl$n_rep <- length(unique(prototype$rep))
psc_only <- prototype[prototype$method == "PSC proximity-weighted", ]
acc_rate <- aggregate(accepted ~ overlap + method, data = psc_only, FUN = mean)
names(acc_rate)[names(acc_rate) == "accepted"] <- "acceptance"
acc_n <- aggregate(accepted ~ overlap + method, data = psc_only, FUN = sum)
names(acc_n)[names(acc_n) == "accepted"] <- "n_accepted"
summary_tbl <- merge(summary_tbl, acc_rate, by = c("overlap", "method"), all.x = TRUE)
summary_tbl <- merge(summary_tbl, acc_n, by = c("overlap", "method"), all.x = TRUE)
summary_tbl <- summary_tbl[order(summary_tbl$overlap, summary_tbl$method), ]
rownames(summary_tbl) <- NULL

sum_out <- file.path(sim_dir, "data/prototype_single_arm_summary.csv")
write.csv(prototype, reps_out, row.names = FALSE)
write.csv(summary_tbl, sum_out, row.names = FALSE)
print(summary_tbl[, c(
  "overlap", "method", "mean_est", "median_est", "mc_sd", "mean_bias",
  "coverage", "mean_sum_w", "n_accepted", "acceptance"
)])
cat("Saved", reps_out, "\n")
cat("Saved", sum_out, "\n")
