# Internal validation for PSC counterfactual models (CFM).
# Gate uses control-data checks only (not treated-arm support / overlap).

#' Landmark Kaplan-Meier survival at time t (right-censored Surv).
sim_km_surv <- function(time, status, t) {
  if (length(time) < 5L) {
    return(NA_real_)
  }
  sf <- survival::survfit(survival::Surv(time, status) ~ 1)
  i <- findInterval(t, sf$time)
  if (i < 1L) {
    return(1)
  }
  as.numeric(sf$surv[i])
}

#' Mean predicted survival from a flexsurvreg / flexsurvspline at time t.
sim_cfm_mean_surv <- function(fit, newdata, t) {
  pr <- summary(fit, newdata = newdata, t = t, type = "survival", tidy = TRUE)
  mean(pr$est, na.rm = TRUE)
}

#' Harrell C-index using linear predictor from a flexsurv model.
sim_cfm_cindex <- function(fit, newdata, time, status) {
  pr <- tryCatch(
    predict(fit, newdata = newdata, type = "lp"),
    error = function(e) NULL
  )
  if (is.null(pr)) {
    return(NA_real_)
  }
  if (is.data.frame(pr) || inherits(pr, "tbl_df")) {
    col <- intersect(c(".pred_link", ".pred", "pred", "estimate"), names(pr))
    if (!length(col)) {
      return(NA_real_)
    }
    lp <- as.numeric(pr[[col[1]]])
  } else {
    lp <- as.numeric(pr)
  }
  if (!length(lp) || any(!is.finite(lp))) {
    return(NA_real_)
  }
  # flexsurv predict(type="lp") returns a link that is inverted vs Cox risk
  # scores for concordance; use -lp so higher risk → shorter survival.
  risk <- -lp
  suppressWarnings(
    as.numeric(survival::concordance(survival::Surv(time, status) ~ risk)$concordance)
  )
}

#' k-fold internal validation of a CFM formula on control data.
#'
#' Metrics (averaged over folds):
#'   cal_abs_t*  — |mean predicted S(t) − KM(t)| on held-out fold
#'   cindex      — discrimination (concordance) on held-out fold
#'
#' Pass if max calibration absolute error ≤ cal_tol and mean C ≥ cindex_min.
sim_cfm_internal_validate <- function(
    contDat,
    formula = s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5,
    k = 5L,
    landmarks = c(2, 4),
    cal_tol = 0.05,
    cindex_min = 0.55,
    seed = 1L) {

  n <- nrow(contDat)
  if (n < max(30L, 2L * k)) {
    return(list(
      pass = FALSE,
      reason = "too_few_controls",
      cal_abs = NA_real_,
      cindex = NA_real_,
      cal_by_t = setNames(rep(NA_real_, length(landmarks)), paste0("t", landmarks)),
      n = n
    ))
  }

  set.seed(seed)
  folds <- sample(rep_len(seq_len(k), n))
  cal_mat <- matrix(NA_real_, nrow = k, ncol = length(landmarks))
  colnames(cal_mat) <- paste0("t", landmarks)
  c_fold <- rep(NA_real_, k)

  for (f in seq_len(k)) {
    train <- contDat[folds != f, , drop = FALSE]
    test <- contDat[folds == f, , drop = FALSE]
    if (sum(train$status) < 5L || nrow(test) < 5L) {
      next
    }
    fit <- tryCatch(
      flexsurv::flexsurvspline(formula, data = train, k = 3),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      next
    }
    for (j in seq_along(landmarks)) {
      t <- landmarks[j]
      km <- sim_km_surv(test$time, test$status, t)
      pr <- tryCatch(
        sim_cfm_mean_surv(fit, test, t),
        error = function(e) NA_real_
      )
      if (is.finite(km) && is.finite(pr)) {
        cal_mat[f, j] <- abs(pr - km)
      }
    }
    c_fold[f] <- tryCatch(
      sim_cfm_cindex(fit, test, test$time, test$status),
      error = function(e) NA_real_
    )
  }

  cal_by_t <- colMeans(cal_mat, na.rm = TRUE)
  cal_abs <- max(cal_by_t, na.rm = TRUE)
  if (!is.finite(cal_abs)) {
    cal_abs <- Inf
  }
  cindex <- mean(c_fold, na.rm = TRUE)
  if (!is.finite(cindex)) {
    cindex <- NA_real_
  }

  pass_cal <- is.finite(cal_abs) && cal_abs <= cal_tol
  pass_c <- is.finite(cindex) && cindex >= cindex_min
  pass <- isTRUE(pass_cal && pass_c)
  reason <- if (pass) {
    "ok"
  } else if (!pass_cal && !pass_c) {
    "cal_and_cindex_fail"
  } else if (!pass_cal) {
    "calibration_fail"
  } else {
    "cindex_fail"
  }

  list(
    pass = pass,
    reason = reason,
    cal_abs = cal_abs,
    cindex = cindex,
    cal_by_t = cal_by_t,
    n = n,
    cal_tol = cal_tol,
    cindex_min = cindex_min
  )
}
