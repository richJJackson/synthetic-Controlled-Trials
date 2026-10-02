# Personalised synthetic control with adjusted joint proximity weights.
# A Royston-Parmar counterfactual model in X_1-X_6, with 1 internal knot,
# is fitted to the historical controls and applied to the
# trial patients. Joint Mahalanobis weights are inflated by min(1, w/0.5).

load_psc_for_prototype <- function() {
  if (isTRUE(get0(".psc_prototype_ready", inherits = TRUE))) {
    return(invisible(TRUE))
  }
  suppressPackageStartupMessages({
    library(survival)
    library(flexsurv)
    library(psc)
    library(posterior)
  })
  rnew <- "/Users/richardjackson/Documents/GitHub/pscDevelop/psc dev/psc/Rnew"
  bridge <- new.env(parent = asNamespace("psc"))
  for (nm in getNamespaceExports("posterior")) {
    assign(nm, getExportedValue("posterior", nm), envir = bridge)
  }
  dev_env <- new.env(parent = bridge)
  sys.source(file.path(rnew, "covRef.R"), dev_env)
  sys.source(file.path(rnew, "pscSupportWeights.R"), dev_env)
  sys.source(file.path(rnew, "pscData_subFn.R"), dev_env)
  sys.source(file.path(rnew, "pscData.R"), dev_env)
  sys.source(file.path(rnew, "lik.flexsurvreg.R"), dev_env)
  sys.source(file.path(rnew, "init.R"), dev_env)
  sys.source(file.path(rnew, "pscEst_proposal.R"), dev_env)
  sys.source(file.path(rnew, "pscEst_start.R"), dev_env)
  sys.source(file.path(rnew, "pscEst_subFn.R"), dev_env)
  sys.source(file.path(rnew, "pscEst_update.R"), dev_env)
  sys.source(file.path(rnew, "pscEst.R"), dev_env)
  sys.source(file.path(rnew, "postSummary.R"), dev_env)
  sys.source(file.path(rnew, "pscfit.R"), dev_env)
  eval(quote(
    as_pscCFM_dev <- function(cfm) {
      me <- modelExtract.flexsurvreg(cfm)
      me$cov_ref <- covRef_build_extended(cfm, cov_co = me$cov_co)
      me$datavis <- NULL
      class(me) <- "pscCFM"
      me
    }
  ), envir = dev_env)
  assignInNamespace("pscData", dev_env$pscData, "psc")
  assignInNamespace("pscfit", dev_env$pscfit, "psc")
  assignInNamespace("lik.flexsurvreg", dev_env$lik.flexsurvreg, "psc")
  assign("as_pscCFM_dev", dev_env$as_pscCFM_dev, .GlobalEnv)
  assign("pscSupportWeights", dev_env$pscSupportWeights, .GlobalEnv)
  assign(".psc_prototype_ready", TRUE, .GlobalEnv)
  invisible(TRUE)
}

fit_psc <- function(historical, trial, nsim = 2000L, nchain = 1L, seed = NULL) {
  load_psc_for_prototype()
  if (!is.null(seed)) set.seed(seed)
  historical$s.ob <- survival::Surv(historical$time, historical$status)
  historical$X_4 <- factor(historical$X_4, levels = 1:3)
  historical$X_5 <- factor(historical$X_5, levels = 1:3)
  historical$X_6 <- factor(historical$X_6, levels = 1:3)
  trial$X_4 <- factor(trial$X_4, levels = levels(historical$X_4))
  trial$X_5 <- factor(trial$X_5, levels = levels(historical$X_5))
  trial$X_6 <- factor(trial$X_6, levels = levels(historical$X_6))
  trial$cen <- trial$status
  trial$s.ob <- survival::Surv(trial$time, trial$cen)

  n_hist <- nrow(historical)
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
    nchain = nchain
  )
  cal_beta <- as.numeric(posterior::as_draws_df(cal_fit$draws)$beta_1)
  cal_beta <- cal_beta[is.finite(cal_beta)]
  cal_hr <- if (length(cal_beta) == 0L) NA_real_ else stats::median(exp(cal_beta))
  accepted <- is.finite(cal_hr) && cal_hr >= 0.8 && cal_hr <= 1.2
  if (!accepted) {
    return(data.frame(
      method = "PSC proximity-weighted",
      est = NA_real_,
      se = NA_real_,
      lo = NA_real_,
      hi = NA_real_,
      max_imbalance = NA_real_,
      max_weight = NA_real_,
      sum_w = NA_real_,
      ess = NA_real_,
      weight_on = "trial patients",
      accepted = FALSE,
      calibration_slope = cal_hr,
      stringsAsFactors = FALSE
    ))
  }

  cfm <- flexsurv::flexsurvspline(
    s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5 + X_6,
    data = historical,
    k = 1
  )
  cfmw <- as_pscCFM_dev(cfm)
  pw <- pscSupportWeights(cfmw, trial, method = "joint")
  pw$weights <- pmin(1, pw$weights / 0.5)
  fit <- psc::pscfit(
    cfmw,
    trial,
    nsim = nsim,
    nchain = nchain,
    proximity_weights = pw
  )
  beta <- as.numeric(posterior::as_draws_df(fit$draws)$beta_1)
  beta <- beta[is.finite(beta)]
  q <- stats::quantile(beta, c(0.025, 0.975), names = FALSE)
  w <- pw$weights[is.finite(pw$weights)]
  sum_w <- sum(w)
  data.frame(
    method = "PSC proximity-weighted",
    est = stats::median(beta),
    se = stats::sd(beta),
    lo = q[1],
    hi = q[2],
    max_imbalance = NA_real_,
    max_weight = max(w),
    sum_w = sum_w,
    ess = sum_w^2 / sum(w^2),
    weight_on = "trial patients",
    accepted = TRUE,
    calibration_slope = cal_hr,
    stringsAsFactors = FALSE
  )
}
