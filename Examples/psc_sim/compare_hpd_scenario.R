# Run HPD-interval rerun and compare with the original Wald-based result file.
#
# Usage (from psc_sim/):
#   source("compare_hpd_scenario.R")
#   compare_hpd_scenario(nsim = 200L)

compare_hpd_scenario <- function(
    overlap = "no_overlap",
    sample_size = "small",
    trial_type = "single_arm",
    nsim = 200L,
    result_suffix = "_mcmcfix",
    force_restart = TRUE,
    paths = NULL,
    run = TRUE) {

  if (is.null(paths)) {
    paths <- sim_paths()
  }

  id <- scenario_id(overlap, sample_size, trial_type)
  run_id <- paste0(id, result_suffix)
  orig_file <- sim_scenario_result_path(id, paths)
  new_file <- file.path(paths$results_dir, paste0(run_id, ".RData"))

  if (!file.exists(orig_file)) {
    stop("Original result not found: ", orig_file)
  }

  if (run) {
    if (!exists("sim_run_scenario")) {
      source(file.path(paths$sim_dir, "run_scenario.R"), local = FALSE)
    }
    sim_run_scenario(
      overlap = overlap,
      sample_size = sample_size,
      trial_type = trial_type,
      nsim = nsim,
      force_restart = force_restart,
      resume = !force_restart,
      verbose = TRUE,
      result_suffix = result_suffix,
      update_progress = FALSE
    )
  }

  if (!file.exists(new_file)) {
    stop("New result not found: ", new_file)
  }

  orig_env <- new.env()
  load(orig_file, envir = orig_env)
  new_env <- new.env()
  load(new_file, envir = new_env)

  method_names <- if (trial_type == "single_arm") SA_METHOD_NAMES else HYBRID_METHOD_NAMES
  compare_methods <- c(
    "Pooled (full adj)", "Pooled (partial adj)",
    "PSC full", "PSC partial",
    "PSC full (pw)", "PSC partial (pw)",
    "SC full"
  )
  if (trial_type == "hybrid") {
    compare_methods <- c(
      "Pooled (adj)", "PSC full (combined)", "PSC partial (combined)",
      "PSC full (pw, combined)", "SC full"
    )
  }
  compare_methods <- intersect(compare_methods, method_names)

  orig_res <- orig_env$result$res.array
  new_res <- new_env$result$res.array
  match_n <- min(dim(orig_res)[3], dim(new_res)[3])
  if (match_n < dim(orig_res)[3]) {
    orig_res <- orig_res[, , seq_len(match_n), drop = FALSE]
  }

  orig_sum <- summarise_res_array(
    orig_res,
    orig_env$result$true_beta,
    method_names = method_names
  )
  new_sum <- summarise_res_array(
    new_res,
    new_env$result$true_beta,
    method_names = method_names
  )

  orig_sum$interval <- paste0("Wald (original, first ", match_n, " reps)")
  new_sum$interval <- paste0("PSC rerun (", result_suffix, ", n = ", match_n, ")")
  orig_sum$run_n <- orig_sum$n
  new_sum$run_n <- new_sum$n

  orig_cmp <- orig_sum[orig_sum$method %in% compare_methods, ]
  new_cmp <- new_sum[new_sum$method %in% compare_methods, ]

  side <- merge(
    orig_cmp[, c("method", "run_n", "mean_bias", "acil", "coverage")],
    new_cmp[, c("method", "run_n", "mean_bias", "acil", "coverage")],
    by = "method",
    suffixes = c("_wald", "_new")
  )
  side$delta_coverage <- side$coverage_new - side$coverage_wald
  side$delta_acil <- side$acil_new - side$acil_wald

  out_csv <- file.path(paths$results_dir, paste0(run_id, "_vs_original.csv"))
  write.csv(side, out_csv, row.names = FALSE)

  cat("\nComparison:", id, " (", result_suffix, " vs original)\n", sep = "")
  cat("Original:", orig_file, " (n = ", orig_cmp$run_n[1], ")\n", sep = "")
  cat("Rerun:  ", new_file, " (n = ", new_cmp$run_n[1], ")\n\n", sep = "")
  side_print <- side
  num_cols <- vapply(side_print, is.numeric, logical(1))
  side_print[num_cols] <- lapply(side_print[num_cols], function(x) round(x, 4))
  print(side_print, row.names = FALSE)
  cat("\nSaved:", out_csv, "\n")

  invisible(list(
    original = orig_env$result,
    hpd = new_env$result,
    comparison = side,
    comparison_file = out_csv
  ))
}

if (sys.nframe() == 0L || !interactive()) {
  PSC_SIM_DIR <- Sys.getenv("PSC_SIM_DIR", unset = normalizePath("."))
  assign("PSC_SIM_DIR", PSC_SIM_DIR, envir = .GlobalEnv)
  source(file.path(PSC_SIM_DIR, "setup.R"))
  sim_source_all(PSC_SIM_DIR)
  source(file.path(PSC_SIM_DIR, "run_scenario.R"))
  sim_load_packages()
  compare_hpd_scenario(
    nsim = as.integer(Sys.getenv("HPD_COMPARE_NSIM", "200")),
    result_suffix = Sys.getenv("HPD_COMPARE_SUFFIX", "_mcmcfix"),
    force_restart = as.logical(Sys.getenv("HPD_COMPARE_FORCE_RESTART", "TRUE"))
  )
}
