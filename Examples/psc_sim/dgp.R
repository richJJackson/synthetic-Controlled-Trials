# Data-generating process for PSC simulation super-populations.

sim_lp <- function(X_1, X_2, X_3, X_4, X_5, dgp, include_beta = FALSE) {
  lp <- dgp$gam_0 +
    dgp$gam_1 * X_1 +
    dgp$gam_2 * X_2 +
    dgp$gam_3 * X_3 +
    dgp$gam_4 * X_4 +
    dgp$gam_5 * X_5
  if (include_beta) {
    lp <- lp + dgp$beta
  }
  exp(lp)
}

sim_draw_covariates <- function(n, spec) {
  data.frame(
    X_1 = rnorm(n, spec$X_1["mean"], spec$X_1["sd"]),
    X_2 = rnorm(n, spec$X_2["mean"], spec$X_2["sd"]),
    X_3 = rbinom(n, 1, spec$X_3_p),
    X_4 = rbinom(n, 1, spec$X_4_p),
    X_5 = rbinom(n, 1, spec$X_5_p)
  )
}

sim_apply_censoring_control <- function(dat, dgp) {
  dat$cen <- as.numeric(!(dat$centm < dat$tm))
  dat$tm[which(dat$cen == 0)] <- dat$centm[which(dat$cen == 0)]
  dat$cenc <- as.numeric(!(dat$centm < dat$tmc))
  dat$tmc[which(dat$cenc == 0)] <- dat$centm[which(dat$cenc == 0)]
  dat
}

sim_apply_censoring_trial <- function(dat, dgp) {
  dat$cen <- as.numeric(!(dat$centm < dat$tm))
  dat$tm[which(dat$cen == 0)] <- dat$centm[which(dat$cen == 0)]
  dat$cen <- as.numeric(!(dgp$admin_censor_years < dat$tm))
  dat$tm <- pmin(dgp$admin_censor_years, dat$tm)
  dat$cenc <- as.numeric(!(dat$centm < dat$tmc))
  dat$tmc[which(dat$cenc == 0)] <- dat$centm[which(dat$cenc == 0)]
  dat$cenc <- as.numeric(!(dgp$admin_censor_years < dat$tmc))
  dat$tmc <- pmin(dgp$admin_censor_years, dat$tmc)
  dat
}

#' Generate one super-population with survival outcomes.
#'
#' @param n Population size.
#' @param spec Covariate spec from SIM_OVERLAP_SPECS.
#' @param dgp DGP parameters (SIM_DGP).
#' @param role \code{"control"} or \code{"trial"}.
sim_generate_population <- function(n, spec, dgp, role = c("control", "trial")) {
  role <- match.arg(role)
  cov <- sim_draw_covariates(n, spec)
  lp <- sim_lp(cov$X_1, cov$X_2, cov$X_3, cov$X_4, cov$X_5, dgp, include_beta = role == "trial")
  lpc <- sim_lp(cov$X_1, cov$X_2, cov$X_3, cov$X_4, cov$X_5, dgp, include_beta = FALSE)
  tm <- rexp(n, lp)
  tmc <- rexp(n, lpc)
  centm <- rexp(n, dgp$censor_rate)
  dat <- data.frame(cov, tm = tm, tmc = tmc, centm = centm)
  if (role == "control") {
    sim_apply_censoring_control(dat, dgp)
  } else {
    sim_apply_censoring_trial(dat, dgp)
  }
}

#' Generate and return all super-populations for the simulation study.
sim_generate_all_populations <- function(dgp = SIM_DGP, overlap_specs = SIM_OVERLAP_SPECS, seed = dgp$seed_pop %||% 2026L) {
  set.seed(seed)
  control_spec <- overlap_specs$large_overlap
  populations <- list(
    control = sim_generate_population(dgp$n_pop, control_spec, dgp, role = "control")
  )
  for (nm in names(overlap_specs)) {
    populations[[paste0("trial_", nm)]] <- sim_generate_population(
      dgp$n_pop, overlap_specs[[nm]], dgp, role = "trial"
    )
  }
  attr(populations, "dgp") <- dgp
  attr(populations, "overlap_specs") <- overlap_specs
  attr(populations, "seed") <- seed
  populations
}

trial_population_name <- function(overlap) {
  paste0("trial_", overlap)
}

#' Sample analysis cohorts for one simulation replicate.
sim_sample_cohorts <- function(populations, overlap, n_cont, n_trt, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  control <- populations$control
  trial_name <- trial_population_name(overlap)
  if (!trial_name %in% names(populations)) {
    stop("Missing trial population: ", trial_name)
  }
  trial <- populations[[trial_name]]
  list(
    contDat = control[sample.int(nrow(control), n_cont, replace = FALSE), , drop = FALSE],
    trtDat = trial[sample.int(nrow(trial), n_trt, replace = FALSE), , drop = FALSE]
  )
}
