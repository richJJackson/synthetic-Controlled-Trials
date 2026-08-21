# Run one simulation scenario (single-arm or hybrid).

sim_run_scenario <- function(
    overlap,
    sample_size,
    trial_type,
    populations = NULL,
    populations_file = NULL,
    nsim = NULL,
    defaults = SIM_DEFAULTS,
    paths = sim_paths(),
    save_results = defaults$save_results,
    verbose = defaults$verbose,
    resume = TRUE,
    force_restart = FALSE,
    checkpoint_every = defaults$checkpoint_every,
    result_suffix = NULL,
    update_progress = NULL) {

  overlap <- match.arg(overlap, SIM_OVERLAP_LEVELS)
  sample_size <- match.arg(sample_size, names(SIM_SAMPLE_SIZES))
  trial_type <- match.arg(trial_type, SIM_TRIAL_TYPES)

  if (is.null(populations)) {
    populations_file <- populations_file %||% file.path(paths$data_dir, "populations.RData")
    if (!file.exists(populations_file)) {
      stop("Populations not found. Run generate_populations.R first: ", populations_file)
    }
    pop_env <- new.env()
    load(populations_file, envir = pop_env)
    populations <- pop_env$populations
  }

  sizes <- SIM_SAMPLE_SIZES[[sample_size]]
  n_cont <- sizes$n_cont
  n_trt <- sizes$n_trt
  if (is.null(nsim)) {
    nsim <- if (trial_type == "hybrid") defaults$nsim_hybrid else defaults$nsim
  }

  n_methods <- if (trial_type == "single_arm") length(SA_METHOD_NAMES) else length(HYBRID_METHOD_NAMES)
  id <- scenario_id(overlap, sample_size, trial_type)
  run_id <- if (is.null(result_suffix) || !nzchar(result_suffix)) {
    id
  } else {
    paste0(id, result_suffix)
  }
  if (is.null(update_progress)) {
    update_progress <- identical(run_id, id)
  }
  result_file <- file.path(paths$results_dir, paste0(run_id, ".RData"))

  if (update_progress && resume && !force_restart &&
      sim_scenario_is_complete(id, nsim = nsim, paths = paths)) {
    if (verbose) {
      cat("Skipping completed scenario:", id, "\n")
    }
    env <- new.env()
    load(result_file, envir = env)
    return(invisible(env$result))
  }

  res.array <- array(NA_real_, dim = c(n_methods, SIM_RES_NCOL, nsim))
  sim_fail <- rep(NA_character_, nsim)
  start_replicate <- 1L
  started_at <- sim_now_utc()
  t0 <- proc.time()[["elapsed"]]

  if (resume && !force_restart) {
    state <- if (update_progress) {
      sim_load_resume_state(id, nsim = nsim, paths = paths)
    } else {
      sim_load_alt_result_state(run_id, nsim = nsim, paths = paths)
    }
    if (!is.null(state)) {
      merged <- sim_merge_resume_state(state, n_methods, nsim)
      res.array <- merged$res.array
      sim_fail <- merged$sim_fail
      start_replicate <- merged$next_replicate
      started_at <- merged$started_at %||% started_at
      if (verbose) {
        cat(
          "Resuming", id, "from", merged$source,
          ", replicate ", start_replicate, "/", nsim,
          " (", merged$attempted, " attempted)\n",
          sep = ""
        )
      }
    }
  } else if (force_restart) {
    if (update_progress) {
      sim_clear_checkpoint(id, paths)
    } else {
      sim_clear_checkpoint(run_id, paths)
    }
    if (sim_safe_file_exists(result_file)) {
      file.remove(result_file)
    }
  }

  if (start_replicate > nsim) {
    done_so_far <- sim_count_done_replicates(res.array)
    if (done_so_far == 0L) {
      # Full pass ran but every replicate failed — restart after fixing the bug.
      if (verbose) {
        cat("Restarting failed scenario (0 successes):", id, "\n")
      }
      res.array <- array(NA_real_, dim = c(n_methods, SIM_RES_NCOL, nsim))
      sim_fail <- rep(NA_character_, nsim)
      start_replicate <- 1L
      started_at <- sim_now_utc()
      t0 <- proc.time()[["elapsed"]]
      sim_clear_checkpoint(if (update_progress) id else run_id, paths)
    } else {
      if (verbose) {
        cat("Scenario already complete:", id, "\n")
      }
      env <- new.env()
      load(result_file, envir = env)
      return(invisible(env$result))
    }
  }

  if (update_progress) {
    sim_progress_update(
      id,
      status = "running",
      nsim_done = start_replicate - 1L,
      next_replicate = start_replicate,
      started_at = started_at,
      paths = paths
    )
  }

  if (verbose) {
    cat("Scenario:", run_id, " (nsim =", nsim, ")\n", sep = "")
    cores <- sim_parallel_cores(defaults)
    if (cores > 1L) {
      cat("Parallel replicates:", cores, "cores per batch\n")
    }
  }

  run_one <- if (trial_type == "single_arm") sim_run_single_arm else sim_run_hybrid

  run_one_replicate <- function(s) {
    cohorts <- sim_sample_cohorts(
      populations, overlap, n_cont, n_trt,
      seed = defaults$stan_seed + s
    )
    run_one(cohorts$contDat, cohorts$trtDat, defaults = defaults)
  }

  run_one_replicate_safe <- function(s) {
    tryCatch(
      list(s = s, ok = TRUE, out = run_one_replicate(s)),
      error = function(e) {
        list(s = s, ok = FALSE, msg = conditionMessage(e))
      }
    )
  }

  apply_replicate_results <- function(results) {
    for (res in results) {
      if (res$ok) {
        res.array[, , res$s] <<- res$out
      } else {
        sim_fail[res$s] <<- res$msg
        warning(id, " replicate ", res$s, " failed: ", res$msg, call. = FALSE)
      }
    }
  }

  maybe_checkpoint <- function(last_s) {
    elapsed <- proc.time()[["elapsed"]] - t0
    done <- sim_count_done_replicates(res.array)
    attempted <- sim_count_attempted_replicates(res.array, sim_fail)
    sec_per_rep <- elapsed / max(1L, attempted - (start_replicate - 1L))
    if (save_results && (last_s %% checkpoint_every == 0L || last_s == nsim)) {
      sim_save_checkpoint(
        if (update_progress) id else run_id,
        res.array = res.array,
        sim_fail = sim_fail,
        nsim_done = done,
        next_replicate = attempted + 1L,
        nsim_target = nsim,
        started_at = started_at,
        elapsed_sec = elapsed,
        sec_per_replicate = sec_per_rep,
        paths = paths,
        update_progress = update_progress
      )
      if (verbose && last_s != nsim) {
        cat("  checkpoint saved (", done, "/", nsim, ")\n", sep = "")
      }
    }
    if (verbose && (last_s %% max(1L, nsim %/% 10L) == 0L || last_s == nsim)) {
      cat("  replicate", last_s, "/", nsim, "\n")
    }
  }

  cores <- sim_parallel_cores(defaults)
  use_parallel <- cores > 1L && .Platform$OS.type == "unix"

  if (use_parallel) {
    s_cursor <- start_replicate
    while (s_cursor <= nsim) {
      batch_end <- min(s_cursor + cores - 1L, nsim)
      batch_idx <- s_cursor:batch_end
      batch_results <- parallel::mclapply(
        batch_idx,
        run_one_replicate_safe,
        mc.cores = length(batch_idx)
      )
      apply_replicate_results(batch_results)
      maybe_checkpoint(max(batch_idx))
      s_cursor <- batch_end + 1L
    }
  } else {
    for (s in start_replicate:nsim) {
      res <- run_one_replicate_safe(s)
      apply_replicate_results(list(res))
      maybe_checkpoint(s)
    }
  }

  method_names <- if (trial_type == "single_arm") SA_METHOD_NAMES else HYBRID_METHOD_NAMES
  summary_tab <- summarise_res_array(
    res.array,
    true_val = SIM_DGP$beta,
    method_names = method_names
  )

  result <- list(
    scenario_id = id,
    run_id = run_id,
    overlap = overlap,
    sample_size = sample_size,
    trial_type = trial_type,
    n_cont = n_cont,
    n_trt = n_trt,
    nsim = nsim,
    res.array = res.array,
    sim_fail = sim_fail,
    summary = summary_tab,
    true_beta = SIM_DGP$beta,
    completed_at = sim_now_utc()
  )

  if (save_results) {
    sim_ensure_dirs(paths)
    save(result, file = result_file)
    attempted <- sim_count_attempted_replicates(res.array, sim_fail)
    if (attempted >= nsim) {
      sim_clear_checkpoint(if (update_progress) id else run_id, paths)
    }
    if (update_progress) {
      sim_progress_update(
        id,
        status = if (attempted >= nsim) {
          if (sim_count_done_replicates(res.array) == 0L) "failed" else "completed"
        } else {
          "partial"
        },
        nsim_done = sim_count_done_replicates(res.array),
        n_failed = sum(!is.na(sim_fail)),
        next_replicate = attempted + 1L,
        completed_at = if (attempted >= nsim && sim_count_done_replicates(res.array) > 0L) {
          result$completed_at
        } else {
          NULL
        },
        elapsed_sec = proc.time()[["elapsed"]] - t0,
        sec_per_replicate = (proc.time()[["elapsed"]] - t0) / max(1L, nsim - start_replicate + 1L),
        paths = paths
      )
    }
    if (verbose) {
      cat("Saved:", result_file, "\n")
    }
    result$output_file <- result_file
  }

  invisible(result)
}

sim_run_batch <- function(
    grid = scenario_grid(),
    populations_file = NULL,
    resume = TRUE,
    force_restart = FALSE,
    defaults = SIM_DEFAULTS,
    paths = sim_paths(),
    verbose = TRUE) {

  sim_load_manifest(paths, create = TRUE)
  populations_file <- populations_file %||% file.path(paths$data_dir, "populations.RData")

  for (i in seq_len(nrow(grid))) {
    row <- grid[i, ]
    id <- scenario_id(row$overlap, row$sample_size, row$trial_type)
    if (verbose) {
      cat("\n[", i, "/", nrow(grid), "] ", id, "\n", sep = "")
    }
    sim_run_scenario(
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
    sim_combine_summaries(paths)
  }

  sim_print_progress(paths)
  invisible(sim_sync_manifest_from_disk(paths = paths))
}
