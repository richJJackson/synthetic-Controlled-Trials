# Patch PSC method columns in existing scenario results (corrected MCMC only).
#
# Usage (from psc_sim/):
#   source("patch_psc_scenarios.R")
#   patch_psc_scenario("no_overlap", "small", "single_arm")
#   patch_psc_all_scenarios()
#   Rscript patch_psc_scenarios.R

sim_psc_patch_checkpoint_path <- function(id, paths = sim_paths()) {
  file.path(paths$results_dir, "checkpoints", paste0(id, "_psc_patch_checkpoint.rds"))
}

sim_backup_pre_patch <- function(result_file, id, paths = sim_paths()) {
  backup <- file.path(paths$results_dir, paste0(id, "_pre_psc_patch.RData"))
  if (!file.exists(backup) && file.exists(result_file)) {
    file.copy(result_file, backup)
    cat("Backup saved:", backup, "\n")
  }
  invisible(backup)
}

patch_psc_scenario <- function(
    overlap,
    sample_size,
    trial_type,
    populations = NULL,
    populations_file = NULL,
    defaults = SIM_DEFAULTS,
    paths = sim_paths(),
    verbose = TRUE,
    resume = TRUE,
    force_restart = FALSE,
    checkpoint_every = defaults$checkpoint_every,
    backup = TRUE) {

  overlap <- match.arg(overlap, SIM_OVERLAP_LEVELS)
  sample_size <- match.arg(sample_size, names(SIM_SAMPLE_SIZES))
  trial_type <- match.arg(trial_type, SIM_TRIAL_TYPES)

  if (is.null(populations)) {
    populations_file <- populations_file %||% file.path(paths$data_dir, "populations.RData")
    if (!file.exists(populations_file)) {
      stop("Populations not found: ", populations_file)
    }
    pop_env <- new.env()
    load(populations_file, envir = pop_env)
    populations <- pop_env$populations
  }

  sizes <- SIM_SAMPLE_SIZES[[sample_size]]
  id <- scenario_id(overlap, sample_size, trial_type)
  result_file <- sim_scenario_result_path(id, paths)
  ck_file <- sim_psc_patch_checkpoint_path(id, paths)

  if (!file.exists(result_file)) {
    if (verbose) {
      cat("Skipping (no result file):", id, "\n")
    }
    return(invisible(NULL))
  }

  if (backup) {
    sim_backup_pre_patch(result_file, id, paths)
  }

  env <- new.env()
  load(result_file, envir = env)
  result <- env$result
  nsim <- result$nsim
  res.array <- sim_expand_res_array(result$res.array)
  sim_fail <- result$sim_fail
  psc_rows <- sim_psc_method_rows(trial_type)
  method_names <- if (trial_type == "single_arm") SA_METHOD_NAMES else HYBRID_METHOD_NAMES

  start_s <- 1L
  if (resume && !force_restart && file.exists(ck_file)) {
    ck <- readRDS(ck_file)
    if (!is.null(ck$res.array)) {
      res.array <- sim_expand_res_array(ck$res.array, nsim = nsim)
    }
    if (!is.null(ck$sim_fail)) {
      sim_fail <- ck$sim_fail
    }
    start_s <- ck$next_replicate %||% 1L
    if (verbose) {
      cat("Resuming PSC patch:", id, "from replicate", start_s, "/", nsim, "\n")
    }
  } else if (force_restart && file.exists(ck_file)) {
    file.remove(ck_file)
  }

  if (start_s > nsim) {
    if (verbose) {
      cat("PSC patch already complete:", id, "\n")
    }
    return(invisible(result))
  }

  run_psc <- if (trial_type == "single_arm") sim_run_psc_single_arm else sim_run_psc_hybrid

  run_one <- function(s) {
    cohorts <- sim_sample_cohorts(
      populations, overlap, sizes$n_cont, sizes$n_trt,
      seed = defaults$stan_seed + s
    )
    run_psc(cohorts$contDat, cohorts$trtDat, defaults = defaults)
  }

  run_one_safe <- function(s) {
    tryCatch(
      list(s = s, ok = TRUE, out = run_one(s)),
      error = function(e) list(s = s, ok = FALSE, msg = conditionMessage(e))
    )
  }

  apply_results <- function(results) {
    for (res in results) {
      if (res$ok) {
        res.array[psc_rows, , res$s] <<- res$out
      } else {
        warning(id, " PSC patch replicate ", res$s, " failed: ", res$msg, call. = FALSE)
      }
    }
  }

  save_checkpoint <- function(next_s) {
    ck <- list(
      scenario_id = id,
      res.array = res.array,
      sim_fail = sim_fail,
      next_replicate = next_s,
      nsim_target = nsim,
      patched_at = sim_now_utc()
    )
    dir.create(dirname(ck_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(ck, ck_file)
  }

  save_result <- function() {
    result$res.array <- res.array
    result$sim_fail <- sim_fail
    result$summary <- summarise_res_array(
      res.array,
      result$true_beta,
      method_names = method_names
    )
    result$psc_patched_at <- sim_now_utc()
    result$psc_mcmc_version <- "2.2.0"
    save(result, file = result_file)
    write.csv(
      result$summary,
      file.path(paths$results_dir, paste0(id, "_summary.csv")),
      row.names = FALSE
    )
  }

  if (verbose) {
    cat("PSC patch:", id, " (replicates ", start_s, ":", nsim, ")\n", sep = "")
    cores <- sim_parallel_cores(defaults)
    if (cores > 1L) {
      cat("Parallel replicates:", cores, "cores per batch\n")
    }
  }

  t0 <- proc.time()[["elapsed"]]
  cores <- sim_parallel_cores(defaults)
  use_parallel <- cores > 1L && .Platform$OS.type == "unix"

  if (use_parallel) {
    s_cursor <- start_s
    while (s_cursor <= nsim) {
      batch_end <- min(s_cursor + cores - 1L, nsim)
      batch_idx <- s_cursor:batch_end
      batch_results <- parallel::mclapply(
        batch_idx,
        run_one_safe,
        mc.cores = length(batch_idx)
      )
      apply_results(batch_results)
      if (verbose && (batch_end %% max(1L, nsim %/% 10L) == 0L || batch_end == nsim)) {
        cat("  replicate", batch_end, "/", nsim, "\n")
      }
      if (batch_end %% checkpoint_every == 0L || batch_end == nsim) {
        save_checkpoint(batch_end + 1L)
        save_result()
      }
      s_cursor <- batch_end + 1L
    }
  } else {
    for (s in start_s:nsim) {
      apply_results(list(run_one_safe(s)))
      if (verbose && (s %% max(1L, nsim %/% 10L) == 0L || s == nsim)) {
        cat("  replicate", s, "/", nsim, "\n")
      }
      if (s %% checkpoint_every == 0L || s == nsim) {
        save_checkpoint(s + 1L)
        save_result()
      }
    }
  }

  if (file.exists(ck_file)) {
    file.remove(ck_file)
  }
  save_result()
  sim_combine_summaries(paths)

  if (verbose) {
    elapsed <- proc.time()[["elapsed"]] - t0
    cat("Saved patched:", result_file, " (", round(elapsed, 1), "s)\n", sep = "")
  }

  invisible(result)
}

patch_psc_all_scenarios <- function(
    grid = scenario_grid(),
    populations_file = NULL,
    resume = TRUE,
    force_restart = FALSE,
    defaults = SIM_DEFAULTS,
    paths = sim_paths(),
    verbose = TRUE) {

  populations_file <- populations_file %||% file.path(paths$data_dir, "populations.RData")

  for (i in seq_len(nrow(grid))) {
    row <- grid[i, ]
    id <- scenario_id(row$overlap, row$sample_size, row$trial_type)
    if (verbose) {
      cat("\n[", i, "/", nrow(grid), "] ", id, "\n", sep = "")
    }
    if (!file.exists(sim_scenario_result_path(id, paths))) {
      if (verbose) {
        cat("  skip (no result file)\n")
      }
      next
    }
    patch_psc_scenario(
      overlap = row$overlap,
      sample_size = row$sample_size,
      trial_type = row$trial_type,
      populations_file = populations_file,
      defaults = defaults,
      paths = paths,
      verbose = verbose,
      resume = resume,
      force_restart = force_restart
    )
  }

  if (verbose) {
    cat("\nPSC patch complete. Combined summary updated.\n")
  }
  invisible(TRUE)
}

if (sys.nframe() == 0L || !interactive()) {
  PSC_SIM_DIR <- Sys.getenv("PSC_SIM_DIR", unset = normalizePath("."))
  assign("PSC_SIM_DIR", PSC_SIM_DIR, envir = .GlobalEnv)
  source(file.path(PSC_SIM_DIR, "setup.R"))
  sim_source_all(PSC_SIM_DIR)
  sim_load_packages()
  patch_psc_all_scenarios(
    resume = as.logical(Sys.getenv("PSC_PATCH_RESUME", "TRUE")),
    force_restart = as.logical(Sys.getenv("PSC_PATCH_FORCE_RESTART", "FALSE"))
  )
}
