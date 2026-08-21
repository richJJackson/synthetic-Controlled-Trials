# Summarise res.array objects (estimate + SE + interval bounds per replicate).

summarise_res_array <- function(res, true_val, method_names = NULL, z = qnorm(0.975)) {
  res <- sim_expand_res_array(res, n_methods = dim(res)[1], nsim = dim(res)[3])
  nr <- dim(res)[1]
  nsim <- dim(res)[3]
  est <- matrix(res[, 1, ], nrow = nr, ncol = nsim)
  se <- matrix(res[, 2, ], nrow = nr, ncol = nsim)
  lo <- matrix(res[, 3, ], nrow = nr, ncol = nsim)
  hi <- matrix(res[, 4, ], nrow = nr, ncol = nsim)

  if (is.null(method_names)) {
    method_names <- seq_len(nr)
  } else if (length(method_names) != nr) {
    stop("method_names length must match nrow(res).")
  }

  data.frame(
    method = method_names,
    n = apply(!is.na(est), 1, sum),
    mean_est = apply(est, 1, mean, na.rm = TRUE),
    mean_se = apply(se, 1, mean, na.rm = TRUE),
    mean_bias = apply(est, 1, mean, na.rm = TRUE) - true_val,
    acil = apply(hi - lo, 1, mean, na.rm = TRUE),
    coverage = vapply(seq_len(nr), function(i) {
      ok <- !is.na(est[i, ]) & !is.na(lo[i, ]) & !is.na(hi[i, ])
      if (!any(ok)) {
        return(NA_real_)
      }
      mean(true_val >= lo[i, ok] & true_val <= hi[i, ok])
    }, numeric(1)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

print_res_summary <- function(x, title, true_val) {
  cat("\n", title, "\n", paste(rep("-", 60), collapse = ""), "\n", sep = "")
  cat(
    "True log(HR) = ", round(true_val, 4),
    "  (HR = ", round(exp(true_val), 4), ")\n",
    "Coverage uses stored interval limits per replicate ",
    "(95% HPD for PSC; Wald for other methods).\n\n",
    sep = ""
  )
  print(round(x, 4), row.names = FALSE)
  invisible(x)
}
