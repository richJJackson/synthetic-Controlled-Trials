# Master populations for the restarted simulation.
# Three datasets of 50,000. Dataset 1 is the reference.
# Datasets 1 and 2 have a moderate covariate shift.
# Datasets 1 and 3 have a poor covariate shift.
# Continuous covariates are X_1, X_2, X_3.
# Categorical covariates are X_4, X_5, X_6, each with levels 1, 2 and 3.
# All six covariates are drawn from a Gaussian copula. Pairwise correlations
# of the latent variables lie between 0.3 and 0.7. Categorical covariates are
# the latent variables cut at the dataset-specific margins.
# The control hazard depends on these covariates, with the same coefficients
# in every dataset. Datasets 1, 2 and 3 share the intercept, so conditional
# exchangeability given X holds. Datasets 2 and 3 with an outcome shift use
# the same covariate margins and add a positive amount to the log cumulative
# hazard. That shift is smaller for moderate overlap than for poor overlap.

sim_dir <- "~/Documents/GitHub/synthetic-Controlled-Trials/Examples/psc_sim"
sim_dir <- normalizePath(sim_dir, mustWork = TRUE)
data_dir <- file.path(sim_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

n_master <- 50000L
seed_master <- 2026L

# Control cumulative hazard is a Royston-Parmar spline on the log-time
# scale, with one internal knot, and is proportional in the covariates:
#   log H(t) = gam_0 + x b + gamma1 * log(t) + gamma2 * v(log t)
# gamma2 bends the log cumulative hazard, so a Weibull model is short of
# the truth and a one-knot flexsurvspline is the matching family.
# The reference-covariate median is 3.2 years. Treated cumulative hazard
# multiplies by exp(beta). Potential outcomes are censored at the minimum
# of an exponential censoring time and an administrative horizon.
spline_knots <- log(c(0.15, 1.5, 8))
spline_gamma1 <- 1.15
spline_gamma2 <- -0.05
median_ref <- 3.2
gam_0 <- uniroot(function(g0) {
  flexsurv::qsurvspline(
    0.5,
    gamma = c(g0, spline_gamma1, spline_gamma2),
    knots = spline_knots,
    scale = "hazard"
  ) - median_ref
}, c(-4, 0))$root
dgp <- list(
  knots = spline_knots,
  gamma1 = spline_gamma1,
  gamma2 = spline_gamma2,
  gam_0 = gam_0,
  b_X = c(X_1 = 0.12, X_2 = -0.08, X_3 = 0.10),
  b_Xcat = list(
    X_4 = c("2" = 0.10, "3" = 0.25),
    X_5 = c("2" = 0.08, "3" = 0.20),
    X_6 = c("2" = 0.12, "3" = 0.28)
  ),
  beta = log(0.7),
  censor_rate = 0.02,
  admin_years = 6
)

# Latent correlation of (X_1, X_2, X_3, X_4, X_5, X_6).
# Continuous pairs are the observed correlations. Pairs that involve a
# categorical covariate are attenuated once the latent variable is cut.
latent_cor <- matrix(
  c(
    1.00, 0.55, 0.40, 0.70, 0.60, 0.50,
    0.55, 1.00, 0.65, 0.45, 0.70, 0.55,
    0.40, 0.65, 1.00, 0.50, 0.45, 0.70,
    0.70, 0.45, 0.50, 1.00, 0.60, 0.55,
    0.60, 0.70, 0.45, 0.60, 1.00, 0.50,
    0.50, 0.55, 0.70, 0.55, 0.50, 1.00
  ),
  nrow = 6,
  byrow = TRUE
)
latent_chol <- chol(latent_cor)

# Dataset 1 is the reference margin.
# Dataset 2 shifts the margins enough for moderate overlap with dataset 1.
# Dataset 3 shifts them enough for poor overlap with dataset 1.
# outcome_shift is added to the log cumulative hazard of both potential
# outcomes. Zero leaves the hazard unchanged. The shifted datasets repeat
# the moderate and poor margins with a further increase in the hazard.
dataset_specs <- list(
  dataset_1 = list(
    label = "reference",
    overlap_with_1 = "reference",
    outcome_shift = 0,
    X_mean = c(X_1 = 0, X_2 = 0, X_3 = 0),
    X_sd = c(X_1 = 1, X_2 = 1, X_3 = 1),
    Xcat_prob = list(
      X_4 = c(0.55, 0.30, 0.15),
      X_5 = c(0.50, 0.35, 0.15),
      X_6 = c(0.45, 0.35, 0.20)
    )
  ),
  dataset_2 = list(
    label = "moderate_overlap",
    overlap_with_1 = "moderate",
    outcome_shift = 0,
    X_mean = c(X_1 = 1.25, X_2 = 1.05, X_3 = 0.90),
    X_sd = c(X_1 = 1, X_2 = 1, X_3 = 1),
    Xcat_prob = list(
      X_4 = c(0.35, 0.40, 0.25),
      X_5 = c(0.30, 0.40, 0.30),
      X_6 = c(0.30, 0.40, 0.30)
    )
  ),
  dataset_3 = list(
    label = "poor_overlap",
    overlap_with_1 = "poor",
    outcome_shift = 0,
    X_mean = c(X_1 = 2.80, X_2 = 2.40, X_3 = 2.10),
    X_sd = c(X_1 = 1, X_2 = 1, X_3 = 1),
    Xcat_prob = list(
      X_4 = c(0.10, 0.25, 0.65),
      X_5 = c(0.10, 0.30, 0.60),
      X_6 = c(0.15, 0.25, 0.60)
    )
  ),
  dataset_2_shift = list(
    label = "moderate_overlap_outcome_shift",
    overlap_with_1 = "moderate_outcome_shift",
    outcome_shift = 0.20,
    X_mean = c(X_1 = 1.25, X_2 = 1.05, X_3 = 0.90),
    X_sd = c(X_1 = 1, X_2 = 1, X_3 = 1),
    Xcat_prob = list(
      X_4 = c(0.35, 0.40, 0.25),
      X_5 = c(0.30, 0.40, 0.30),
      X_6 = c(0.30, 0.40, 0.30)
    )
  ),
  dataset_3_shift = list(
    label = "poor_overlap_outcome_shift",
    overlap_with_1 = "poor_outcome_shift",
    outcome_shift = 0.45,
    X_mean = c(X_1 = 2.80, X_2 = 2.40, X_3 = 2.10),
    X_sd = c(X_1 = 1, X_2 = 1, X_3 = 1),
    Xcat_prob = list(
      X_4 = c(0.10, 0.25, 0.65),
      X_5 = c(0.10, 0.30, 0.60),
      X_6 = c(0.15, 0.25, 0.60)
    )
  )
)

cut_latent <- function(z, probs) {
  cuts <- qnorm(cumsum(probs[-length(probs)]))
  out <- rep(length(probs), length(z))
  out[z <= cuts[1]] <- 1L
  if (length(cuts) > 1L) {
    for (k in seq(2L, length(cuts))) {
      out[z > cuts[k - 1L] & z <= cuts[k]] <- k
    }
  }
  out
}

draw_dataset <- function(n, spec, dgp, dataset_id) {
  z <- matrix(rnorm(n * 6L), n, 6L) %*% latent_chol
  X_1 <- spec$X_mean["X_1"] + spec$X_sd["X_1"] * z[, 1]
  X_2 <- spec$X_mean["X_2"] + spec$X_sd["X_2"] * z[, 2]
  X_3 <- spec$X_mean["X_3"] + spec$X_sd["X_3"] * z[, 3]
  X_4 <- cut_latent(z[, 4], spec$Xcat_prob$X_4)
  X_5 <- cut_latent(z[, 5], spec$Xcat_prob$X_5)
  X_6 <- cut_latent(z[, 6], spec$Xcat_prob$X_6)

  lp0 <- dgp$gam_0 + spec$outcome_shift +
    dgp$b_X["X_1"] * X_1 +
    dgp$b_X["X_2"] * X_2 +
    dgp$b_X["X_3"] * X_3 +
    dgp$b_Xcat$X_4["2"] * (X_4 == 2L) +
    dgp$b_Xcat$X_4["3"] * (X_4 == 3L) +
    dgp$b_Xcat$X_5["2"] * (X_5 == 2L) +
    dgp$b_Xcat$X_5["3"] * (X_5 == 3L) +
    dgp$b_Xcat$X_6["2"] * (X_6 == 2L) +
    dgp$b_Xcat$X_6["3"] * (X_6 == 3L)

  gamma0 <- cbind(lp0, dgp$gamma1, dgp$gamma2)
  gamma1 <- cbind(lp0 + dgp$beta, dgp$gamma1, dgp$gamma2)
  t0_lat <- flexsurv::qsurvspline(runif(n), gamma = gamma0, knots = dgp$knots, scale = "hazard")
  t1_lat <- flexsurv::qsurvspline(runif(n), gamma = gamma1, knots = dgp$knots, scale = "hazard")
  cens <- pmin(rexp(n, dgp$censor_rate), dgp$admin_years)

  data.frame(
    id = seq_len(n),
    dataset = dataset_id,
    overlap_with_1 = spec$overlap_with_1,
    X_1 = X_1,
    X_2 = X_2,
    X_3 = X_3,
    X_4 = factor(X_4, levels = 1:3),
    X_5 = factor(X_5, levels = 1:3),
    X_6 = factor(X_6, levels = 1:3),
    t0 = pmin(t0_lat, cens),
    event0 = as.integer(t0_lat <= cens),
    t1 = pmin(t1_lat, cens),
    event1 = as.integer(t1_lat <= cens)
  )
}

set.seed(seed_master)
master <- list(
  dataset_1 = draw_dataset(n_master, dataset_specs$dataset_1, dgp, 1L),
  dataset_2 = draw_dataset(n_master, dataset_specs$dataset_2, dgp, 2L),
  dataset_3 = draw_dataset(n_master, dataset_specs$dataset_3, dgp, 3L),
  dataset_2_shift = draw_dataset(n_master, dataset_specs$dataset_2_shift, dgp, 4L),
  dataset_3_shift = draw_dataset(n_master, dataset_specs$dataset_3_shift, dgp, 5L)
)
attr(master, "dgp") <- dgp
attr(master, "dataset_specs") <- dataset_specs
attr(master, "n") <- n_master
attr(master, "seed") <- seed_master

x_names <- c("X_1", "X_2", "X_3")
cat_names <- c("X_4", "X_5", "X_6")
ref <- master$dataset_1
ref_center <- colMeans(ref[, x_names])
ref_cov <- stats::cov(ref[, x_names])
ref_md <- stats::mahalanobis(ref[, x_names], ref_center, ref_cov)
md_cut <- as.numeric(stats::quantile(ref_md, 0.95))

margin_row <- function(dat, name) {
  md <- stats::mahalanobis(dat[, x_names], ref_center, ref_cov)
  row <- data.frame(
    dataset = name,
    n = nrow(dat),
    X_1_mean = mean(dat$X_1),
    X_2_mean = mean(dat$X_2),
    X_3_mean = mean(dat$X_3),
    outcome_shift = dataset_specs[[name]]$outcome_shift,
    joint_overlap_with_1 = mean(md <= md_cut),
    event0_rate = mean(dat$event0),
    median_t0 = median(dat$t0),
    stringsAsFactors = FALSE
  )
  for (cn in cat_names) {
    props <- prop.table(table(dat[[cn]]))
    for (lev in names(props)) {
      row[[paste0(cn, "_p", lev)]] <- as.numeric(props[[lev]])
    }
  }
  row
}

master_summary <- rbind(
  margin_row(master$dataset_1, "dataset_1"),
  margin_row(master$dataset_2, "dataset_2"),
  margin_row(master$dataset_3, "dataset_3"),
  margin_row(master$dataset_2_shift, "dataset_2_shift"),
  margin_row(master$dataset_3_shift, "dataset_3_shift")
)

save(master, file = file.path(data_dir, "master.RData"))
write.csv(
  master_summary,
  file.path(data_dir, "master_summary.csv"),
  row.names = FALSE
)

observed_cor <- function(dat) {
  num <- data.frame(
    X_1 = dat$X_1,
    X_2 = dat$X_2,
    X_3 = dat$X_3,
    X_4 = as.integer(dat$X_4),
    X_5 = as.integer(dat$X_5),
    X_6 = as.integer(dat$X_6)
  )
  stats::cor(num)
}

cat("Saved", file.path(data_dir, "master.RData"), "\n\n")
cat("Control survival and overlap with dataset 1\n")
print(master_summary[, c(
  "dataset", "outcome_shift", "median_t0", "event0_rate", "joint_overlap_with_1"
)])
cat("\nProportions of X_4, X_5, X_6\n")
print(master_summary[, c(
  "dataset",
  "X_4_p1", "X_4_p2", "X_4_p3",
  "X_5_p1", "X_5_p2", "X_5_p3",
  "X_6_p1", "X_6_p2", "X_6_p3"
)])
cat("\nSpline baseline gam_0", dgp$gam_0, "gamma1", dgp$gamma1, "gamma2", dgp$gamma2, "\n")
set.seed(11)
chk <- master$dataset_1[sample.int(nrow(master$dataset_1), 2000L), ]
chk$s.ob <- survival::Surv(chk$t0, chk$event0)
fit0 <- flexsurv::flexsurvspline(
  s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6, data = chk, k = 0)
fit1 <- flexsurv::flexsurvspline(
  s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6, data = chk, k = 1)
cat("Weibull loglik", logLik(fit0), " one-knot loglik", logLik(fit1), "\n")
cat("One-knot baseline coefficients\n")
print(round(coef(fit1), 3))

for (nm in names(master)) {
  cr <- observed_cor(master[[nm]])
  off <- cr[lower.tri(cr)]
  cat(
    "\n", nm, " observed correlation range: ",
    round(min(off), 3), " to ", round(max(off), 3), "\n",
    sep = ""
  )
  print(round(cr, 3))
}
