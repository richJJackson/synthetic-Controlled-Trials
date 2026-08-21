# Gated CFM validation + doubled MCMC (nsim = 10000) experiment.
#
# Fits CFM → internal validation on controls → PSC only if gate passes.
# Stores results separately under results/gated/ for comparison.
#
# Usage (from psc_sim/):
#   Rscript run_gated_cfm_scenario.R
#   source("run_gated_cfm_scenario.R"); run_gated_cfm_scenario()

run_gated_cfm_scenario <- function(
    overlap = "moderate_overlap",
    sample_size = "medium",
    trial_type = "single_arm",
    nsim = 200L,
    psc_nsim = 10000L,
    cal_tol = 0.05,
    cindex_min = 0.55,
    force_restart = TRUE,
    paths = NULL,
    verbose = TRUE) {

  if (is.null(paths)) {
    paths <- sim_paths()
  }
  if (!exists("sim_run_single_arm_gated", mode = "function")) {
    source(file.path(paths$sim_dir, "cfm_validate.R"), local = FALSE)
    source(file.path(paths$sim_dir, "methods_gated.R"), local = FALSE)
  }

  id <- scenario_id(overlap, sample_size, trial_type)
  out_dir <- file.path(paths$results_dir, "gated")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  run_id <- paste0(id, "_gated_mc", psc_nsim)
  result_file <- file.path(out_dir, paste0(run_id, ".RData"))
  ck_file <- file.path(out_dir, paste0(run_id, "_checkpoint.rds"))
  compare_file <- file.path(out_dir, paste0(run_id, "_vs_current.csv"))

  pop_file <- file.path(paths$data_dir, "populations.RData")
  if (!file.exists(pop_file)) {
    stop("Populations not found: ", pop_file)
  }
  pop_env <- new.env()
  load(pop_file, envir = pop_env)
  populations <- pop_env$populations

  sizes <- SIM_SAMPLE_SIZES[[sample_size]]
  defaults <- modifyList(SIM_DEFAULTS, list(psc_nsim = as.integer(psc_nsim)))
  n_methods <- length(GATED_SA_METHOD_NAMES)

  start_s <- 1L
  res.array <- array(NA_real_, dim = c(n_methods, SIM_RES_NCOL, nsim))
  gate_full_pass <- rep(NA, nsim)
  gate_part_pass <- rep(NA, nsim)
  gate_full_reason <- rep(NA_character_, nsim)
  gate_part_reason <- rep(NA_character_, nsim)
  gate_full_cal_abs <- rep(NA_real_, nsim)
  gate_part_cal_abs <- rep(NA_real_, nsim)
  gate_full_cindex <- rep(NA_real_, nsim)
  gate_part_cindex <- rep(NA_real_, nsim)
  sim_fail <- rep(NA_character_, nsim)

  if (!force_restart && file.exists(ck_file)) {
    ck <- readRDS(ck_file)
    res.array <- ck$res.array
    gate_full_pass <- ck$gate_full_pass
    gate_part_pass <- ck$gate_part_pass
    gate_full_reason <- ck$gate_full_reason
    gate_part_reason <- ck$gate_part_reason
    gate_full_cal_abs <- ck$gate_full_cal_abs
    gate_part_cal_abs <- ck$gate_part_cal_abs
    gate_full_cindex <- ck$gate_full_cindex
    gate_part_cindex <- ck$gate_part_cindex
    sim_fail <- ck$sim_fail
    start_s <- ck$next_replicate
    if (verbose) {
      cat("Resuming gated run from replicate", start_s, "/", nsim, "\n")
    }
  }

  if (start_s > nsim) {
    if (verbose) cat("Gated run already complete:", result_file, "\n")
    load(result_file)
    return(invisible(result))
  }

  if (verbose) {
    cat(
      "Gated CFM + MCMC experiment\n",
      "  scenario: ", id, "\n",
      "  nsim: ", nsim, "  psc_nsim: ", psc_nsim, "\n",
      "  cal_tol: ", cal_tol, "  cindex_min: ", cindex_min, "\n",
      "  output: ", result_file, "\n",
      sep = ""
    )
  }

  t0 <- proc.time()[["elapsed"]]
  cores <- sim_parallel_cores(defaults)
  if (verbose && cores > 1L) {
    cat("Parallel replicates:", cores, "cores per batch\n")
  }

  run_one <- function(s) {
    cohorts <- sim_sample_cohorts(
      populations, overlap, sizes$n_cont, sizes$n_trt,
      seed = defaults$stan_seed + s
    )
    sim_run_single_arm_gated(
      cohorts$contDat,
      cohorts$trtDat,
      defaults = defaults,
      cal_tol = cal_tol,
      cindex_min = cindex_min,
      validate_seed = defaults$stan_seed + 1000L + s
    )
  }

  run_one_safe <- function(s) {
    tryCatch(
      list(s = s, ok = TRUE, out = run_one(s)),
      error = function(e) list(s = s, ok = FALSE, msg = conditionMessage(e))
    )
  }

  apply_one <- function(res) {
    s <- res$s
    if (!res$ok) {
      sim_fail[s] <<- res$msg
      warning(run_id, " replicate ", s, " failed: ", res$msg, call. = FALSE)
      return(invisible(NULL))
    }
    out <- res$out
    res.array[, , s] <<- out$res
    gate_full_pass[s] <<- out$gate_full_pass
    gate_part_pass[s] <<- out$gate_part_pass
    gate_full_reason[s] <<- out$gate_full_reason
    gate_part_reason[s] <<- out$gate_part_reason
    gate_full_cal_abs[s] <<- out$gate_full_cal_abs
    gate_part_cal_abs[s] <<- out$gate_part_cal_abs
    gate_full_cindex[s] <<- out$gate_full_cindex
    gate_part_cindex[s] <<- out$gate_part_cindex
  }

  save_ck <- function(next_s) {
    saveRDS(
      list(
        res.array = res.array,
        gate_full_pass = gate_full_pass,
        gate_part_pass = gate_part_pass,
        gate_full_reason = gate_full_reason,
        gate_part_reason = gate_part_reason,
        gate_full_cal_abs = gate_full_cal_abs,
        gate_part_cal_abs = gate_part_cal_abs,
        gate_full_cindex = gate_full_cindex,
        gate_part_cindex = gate_part_cindex,
        sim_fail = sim_fail,
        next_replicate = next_s,
        nsim_target = nsim
      ),
      ck_file
    )
  }

  use_parallel <- cores > 1L && .Platform$OS.type == "unix"
  checkpoint_every <- defaults$checkpoint_every

  if (use_parallel) {
    s_cursor <- start_s
    while (s_cursor <= nsim) {
      batch_end <- min(s_cursor + cores - 1L, nsim)
      batch <- parallel::mclapply(
        s_cursor:batch_end, run_one_safe, mc.cores = batch_end - s_cursor + 1L
      )
      for (res in batch) apply_one(res)
      if (batch_end %% checkpoint_every == 0L || batch_end == nsim) {
        save_ck(batch_end + 1L)
      }
      if (verbose && (batch_end %% max(1L, nsim %/% 10L) == 0L || batch_end == nsim)) {
        cat(
          "  replicate ", batch_end, "/", nsim,
          "  gate_full_pass=", sum(gate_full_pass[1:batch_end], na.rm = TRUE), "/", batch_end,
          "\n",
          sep = ""
        )
      }
      s_cursor <- batch_end + 1L
    }
  } else {
    for (s in start_s:nsim) {
      apply_one(run_one_safe(s))
      if (s %% checkpoint_every == 0L || s == nsim) {
        save_ck(s + 1L)
      }
      if (verbose && (s %% max(1L, nsim %/% 10L) == 0L || s == nsim)) {
        cat("  replicate", s, "/", nsim, "\n")
      }
    }
  }

  summary_all <- summarise_res_array(res.array, SIM_DGP$beta, GATED_SA_METHOD_NAMES)
  # Coverage among gated-in replicates for PSC full
  pass_idx <- which(gate_full_pass %in% TRUE)
  if (length(pass_idx)) {
    summary_gated <- summarise_res_array(
      res.array[, , pass_idx, drop = FALSE],
      SIM_DGP$beta,
      GATED_SA_METHOD_NAMES
    )
  } else {
    summary_gated <- summary_all
    summary_gated[, c("n", "mean_est", "mean_se", "mean_bias", "acil", "coverage")] <- NA
  }

  result <- list(
    scenario_id = id,
    run_id = run_id,
    overlap = overlap,
    sample_size = sample_size,
    trial_type = trial_type,
    nsim = nsim,
    psc_nsim = psc_nsim,
    cal_tol = cal_tol,
    cindex_min = cindex_min,
    true_beta = SIM_DGP$beta,
    method_names = GATED_SA_METHOD_NAMES,
    res.array = res.array,
    gate_full_pass = gate_full_pass,
    gate_part_pass = gate_part_pass,
    gate_full_reason = gate_full_reason,
    gate_part_reason = gate_part_reason,
    gate_full_cal_abs = gate_full_cal_abs,
    gate_part_cal_abs = gate_part_cal_abs,
    gate_full_cindex = gate_full_cindex,
    gate_part_cindex = gate_part_cindex,
    sim_fail = sim_fail,
    summary_all = summary_all,
    summary_gated_in = summary_gated,
    gate_full_pass_rate = mean(gate_full_pass, na.rm = TRUE),
    gate_part_pass_rate = mean(gate_part_pass, na.rm = TRUE),
    elapsed_sec = proc.time()[["elapsed"]] - t0
  )

  save(result, file = result_file)
  if (file.exists(ck_file)) {
    file.remove(ck_file)
  }

  # Compare with first nsim replicates of current (patched) scenario file
  current_file <- sim_scenario_result_path(id, paths)
  compare <- NULL
  if (file.exists(current_file)) {
    env <- new.env()
    load(current_file, envir = env)
    cur <- env$result$res.array
    n_cmp <- min(nsim, dim(cur)[3])
    cur_sum <- summarise_res_array(
      sim_expand_res_array(cur[, , seq_len(n_cmp), drop = FALSE]),
      SIM_DGP$beta,
      SA_METHOD_NAMES
    )
    methods_cmp <- c(
      "Pooled (full adj)", "Pooled (partial adj)",
      "PSC full", "PSC partial", "PSC full (pw)", "PSC partial (pw)"
    )
    compare <- do.call(rbind, lapply(methods_cmp, function(m) {
      g_all <- summary_all[summary_all$method == m, , drop = FALSE]
      g_in <- summary_gated[summary_gated$method == m, , drop = FALSE]
      c_row <- cur_sum[cur_sum$method == m, , drop = FALSE]
      data.frame(
        method = m,
        n_current = if (nrow(c_row)) c_row$n else NA,
        coverage_current = if (nrow(c_row)) c_row$coverage else NA,
        acil_current = if (nrow(c_row)) c_row$acil else NA,
        bias_current = if (nrow(c_row)) c_row$mean_bias else NA,
        n_gated_all = if (nrow(g_all)) g_all$n else NA,
        coverage_gated_all = if (nrow(g_all)) g_all$coverage else NA,
        n_gated_in = if (nrow(g_in)) g_in$n else NA,
        coverage_gated_in = if (nrow(g_in)) g_in$coverage else NA,
        acil_gated_in = if (nrow(g_in)) g_in$acil else NA,
        bias_gated_in = if (nrow(g_in)) g_in$mean_bias else NA,
        stringsAsFactors = FALSE
      )
    }))
    write.csv(compare, compare_file, row.names = FALSE)
  }

  if (verbose) {
    cat("\nSaved:", result_file, "\n")
    cat(
      "Gate pass rates — full CFM: ",
      round(100 * result$gate_full_pass_rate, 1), "%; partial: ",
      round(100 * result$gate_part_pass_rate, 1), "%\n",
      sep = ""
    )
    cat("\nSummary (gated-in PSC only):\n")
    print(round(summary_gated, 4), row.names = FALSE)
    if (!is.null(compare)) {
      cat("\nComparison vs current (first ", n_cmp, " replicates):\n", sep = "")
      print(compare, row.names = FALSE)
      cat("\nWrote:", compare_file, "\n")
    }
    cat("Elapsed:", round(result$elapsed_sec, 1), "s\n")
  }

  invisible(result)
}

if (sys.nframe() == 0L || !interactive()) {
  PSC_SIM_DIR <- Sys.getenv("PSC_SIM_DIR", unset = normalizePath("."))
  assign("PSC_SIM_DIR", PSC_SIM_DIR, envir = .GlobalEnv)
  Sys.setenv(PSC_SIM_DIR = PSC_SIM_DIR)
  source(file.path(PSC_SIM_DIR, "setup.R"))
  sim_source_all(PSC_SIM_DIR)
  source(file.path(PSC_SIM_DIR, "cfm_validate.R"), local = FALSE)
  source(file.path(PSC_SIM_DIR, "methods_gated.R"), local = FALSE)
  sim_load_packages()
  run_gated_cfm_scenario(
    nsim = as.integer(Sys.getenv("GATED_NSIM", "200")),
    psc_nsim = as.integer(Sys.getenv("GATED_PSC_NSIM", "10000")),
    force_restart = as.logical(Sys.getenv("GATED_FORCE_RESTART", "TRUE"))
  )
}
