# Progress tracking, checkpoints, and resume for the simulation study.

sim_manifest_path <- function(paths = sim_paths()) {
  file.path(paths$results_dir, "progress_manifest.rds")
}

sim_checkpoint_dir <- function(paths = sim_paths()) {
  file.path(paths$results_dir, "checkpoints")
}

sim_scenario_result_path <- function(id, paths = sim_paths()) {
  file.path(paths$results_dir, paste0(id, ".RData"))
}

sim_checkpoint_path <- function(id, paths = sim_paths()) {
  file.path(sim_checkpoint_dir(paths), paste0(id, "_checkpoint.rds"))
}

sim_now_utc <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "")
}

sim_elapsed_sec <- function(started_at) {
  if (is.null(started_at) || is.na(started_at) || !nzchar(started_at)) {
    return(NA_real_)
  }
  st <- suppressWarnings(as.POSIXct(started_at, tz = ""))
  if (is.na(st)) {
    return(NA_real_)
  }
  as.numeric(difftime(Sys.time(), st, units = "secs"))
}

sim_count_done_replicates <- function(res.array) {
  if (is.null(res.array) || length(res.array) == 0L) {
    return(0L)
  }
  # Keep replicate dimension even when nsim == 1.
  est <- res.array[, 1, , drop = FALSE]
  # est dim: methods x 1 x nsim
  sum(!is.na(est[1, 1, ]))
}

sim_count_attempted_replicates <- function(res.array, sim_fail = NULL) {
  if (is.null(res.array) || length(res.array) == 0L) {
    return(0L)
  }
  est <- res.array[, 1, , drop = FALSE]
  ok <- !is.na(est[1, 1, ])
  if (!is.null(sim_fail) && length(sim_fail) == length(ok)) {
    ok <- ok | !is.na(sim_fail)
  }
  sum(ok)
}

#' Load the furthest progress for a scenario from result file and/or checkpoint.
sim_load_resume_state <- function(id, nsim = NULL, paths = sim_paths()) {
  result_file <- sim_scenario_result_path(id, paths)
  ck <- sim_load_checkpoint(id, paths)
  candidates <- list()

  if (sim_safe_file_exists(result_file)) {
    env <- new.env()
    load(result_file, envir = env)
    result <- env$result
    target <- nsim %||% result$nsim
    attempted <- sim_count_attempted_replicates(result$res.array, result$sim_fail)
    candidates$result <- list(
      res.array = result$res.array,
      sim_fail = result$sim_fail,
      attempted = attempted,
      done = sim_count_done_replicates(result$res.array),
      next_replicate = attempted + 1L,
      nsim_target = target,
      started_at = NULL,
      source = "result"
    )
  }

  if (!is.null(ck)) {
    target <- nsim %||% ck$nsim_target
    attempted <- sim_count_attempted_replicates(ck$res.array, ck$sim_fail)
    candidates$checkpoint <- list(
      res.array = ck$res.array,
      sim_fail = ck$sim_fail,
      attempted = attempted,
      done = sim_count_done_replicates(ck$res.array),
      next_replicate = ck$next_replicate %||% (attempted + 1L),
      nsim_target = target,
      started_at = ck$started_at,
      source = "checkpoint"
    )
  }

  if (length(candidates) == 0L) {
    return(NULL)
  }

  best_name <- names(candidates)[which.max(vapply(
    candidates,
    function(x) x$attempted,
    integer(1)
  ))]
  state <- candidates[[best_name]]

  # If result file is ahead, drop stale checkpoint so resume cannot regress.
  if (best_name == "result" && !is.null(ck) &&
      state$attempted > candidates$checkpoint$attempted) {
    sim_clear_checkpoint(id, paths)
  }

  state
}

#' Copy loaded progress into full-size arrays for the target replicate count.
sim_merge_resume_state <- function(state, n_methods, nsim) {
  if (is.null(state)) {
    return(list(
      res.array = array(NA_real_, dim = c(n_methods, SIM_RES_NCOL, nsim)),
      sim_fail = rep(NA_character_, nsim),
      next_replicate = 1L
    ))
  }

  res.array <- array(NA_real_, dim = c(n_methods, SIM_RES_NCOL, nsim))
  sim_fail <- rep(NA_character_, nsim)

  loaded_res <- state$res.array
  loaded_fail <- state$sim_fail
  n_loaded <- if (is.null(loaded_res)) 0L else dim(loaded_res)[3]
  n_copy <- min(n_loaded, nsim)

  if (n_copy > 0L) {
    loaded_res <- sim_expand_res_array(loaded_res, n_methods = n_methods, nsim = n_loaded)
    res.array[, , seq_len(n_copy)] <- loaded_res[, , seq_len(n_copy), drop = FALSE]
  }
  if (!is.null(loaded_fail) && length(loaded_fail) > 0L) {
    n_fail <- min(length(loaded_fail), nsim)
    sim_fail[seq_len(n_fail)] <- loaded_fail[seq_len(n_fail)]
  }

  attempted <- sim_count_attempted_replicates(res.array, sim_fail)
  next_replicate <- min(max(1L, attempted + 1L), nsim + 1L)
  if (!is.null(state$next_replicate)) {
    next_replicate <- min(as.integer(state$next_replicate), nsim + 1L)
    next_replicate <- max(next_replicate, attempted + 1L)
  }

  list(
    res.array = res.array,
    sim_fail = sim_fail,
    next_replicate = next_replicate,
    attempted = attempted,
    started_at = state$started_at,
    source = state$source
  )
}

#' Load progress for an alternate result file (suffix runs), not the main manifest.
sim_load_alt_result_state <- function(run_id, nsim = NULL, paths = sim_paths()) {
  result_file <- file.path(paths$results_dir, paste0(run_id, ".RData"))
  ck <- sim_load_checkpoint(run_id, paths)
  candidates <- list()

  if (sim_safe_file_exists(result_file)) {
    env <- new.env()
    load(result_file, envir = env)
    result <- env$result
    target <- nsim %||% result$nsim
    attempted <- sim_count_attempted_replicates(result$res.array, result$sim_fail)
    candidates$result <- list(
      res.array = result$res.array,
      sim_fail = result$sim_fail,
      attempted = attempted,
      done = sim_count_done_replicates(result$res.array),
      next_replicate = attempted + 1L,
      nsim_target = target,
      started_at = NULL,
      source = "result"
    )
  }

  if (!is.null(ck)) {
    target <- nsim %||% ck$nsim_target
    attempted <- sim_count_attempted_replicates(ck$res.array, ck$sim_fail)
    candidates$checkpoint <- list(
      res.array = ck$res.array,
      sim_fail = ck$sim_fail,
      attempted = attempted,
      done = sim_count_done_replicates(ck$res.array),
      next_replicate = ck$next_replicate %||% (attempted + 1L),
      nsim_target = target,
      started_at = ck$started_at,
      source = "checkpoint"
    )
  }

  if (length(candidates) == 0L) {
    return(NULL)
  }

  best_name <- names(candidates)[which.max(vapply(
    candidates,
    function(x) x$attempted,
    integer(1)
  ))]
  candidates[[best_name]]
}

sim_valid_file <- function(path) {
  is.character(path) && length(path) == 1L && !is.na(path) && nzchar(path)
}

sim_safe_file_exists <- function(path) {
  sim_valid_file(path) && file.exists(path)
}

sim_repair_scenario_entry <- function(entry, id, paths = sim_paths()) {
  if (is.null(entry$scenario_id) || !nzchar(entry$scenario_id)) {
    entry$scenario_id <- id
  }
  if (!sim_valid_file(entry$result_file)) {
    entry$result_file <- sim_scenario_result_path(id, paths)
  }
  if (!sim_valid_file(entry$checkpoint_file)) {
    entry$checkpoint_file <- sim_checkpoint_path(id, paths)
  }
  if (is.null(entry$nsim_target)) {
    parsed <- tryCatch(parse_scenario_id(id), error = function(e) NULL)
    if (!is.null(parsed)) {
      entry$nsim_target <- if (parsed$trial_type == "hybrid") {
        SIM_DEFAULTS$nsim_hybrid
      } else {
        SIM_DEFAULTS$nsim
      }
    }
  }
  entry
}

sim_init_manifest <- function(grid = scenario_grid(), defaults = SIM_DEFAULTS, paths = sim_paths()) {
  sim_ensure_dirs(paths)
  dir.create(sim_checkpoint_dir(paths), recursive = TRUE, showWarnings = FALSE)

  scenarios <- list()
  for (i in seq_len(nrow(grid))) {
    row <- grid[i, ]
    id <- scenario_id(row$overlap, row$sample_size, row$trial_type)
    nsim <- if (row$trial_type == "hybrid") defaults$nsim_hybrid else defaults$nsim
    scenarios[[id]] <- list(
      scenario_id = id,
      overlap = row$overlap,
      sample_size = row$sample_size,
      trial_type = row$trial_type,
      nsim_target = nsim,
      nsim_done = 0L,
      n_failed = 0L,
      status = "pending",
      started_at = NA_character_,
      updated_at = NA_character_,
      completed_at = NA_character_,
      next_replicate = 1L,
      result_file = sim_scenario_result_path(id, paths),
      checkpoint_file = sim_checkpoint_path(id, paths),
      elapsed_sec = NA_real_,
      sec_per_replicate = NA_real_
    )
  }

  list(
    version = 1L,
    study = "psc_sim",
    created_at = sim_now_utc(),
    updated_at = sim_now_utc(),
    defaults = defaults,
    scenarios = scenarios
  )
}

sim_load_manifest <- function(paths = sim_paths(), create = TRUE) {
  manifest_file <- sim_manifest_path(paths)
  if (!sim_valid_file(manifest_file)) {
    stop("Invalid manifest path: ", manifest_file)
  }
  if (file.exists(manifest_file)) {
    manifest <- readRDS(manifest_file)
    bad_names <- !nzchar(names(manifest$scenarios))
    if (any(bad_names)) {
      warning(
        "Manifest has invalid scenario names. Run: source('repair_manifest.R')",
        call. = FALSE
      )
      manifest$scenarios <- manifest$scenarios[!bad_names]
    }
    for (id in names(manifest$scenarios)) {
      manifest$scenarios[[id]] <- sim_repair_scenario_entry(
        manifest$scenarios[[id]], id, paths
      )
    }
    manifest
  } else if (create) {
    manifest <- sim_init_manifest(paths = paths)
    sim_save_manifest(manifest, paths)
    manifest
  } else {
    NULL
  }
}

sim_save_manifest <- function(manifest, paths = sim_paths()) {
  manifest$updated_at <- sim_now_utc()
  saveRDS(manifest, sim_manifest_path(paths))
  invisible(manifest)
}

sim_sync_manifest_from_disk <- function(manifest = NULL, paths = sim_paths()) {
  if (is.null(manifest)) {
    manifest <- sim_load_manifest(paths, create = TRUE)
  }

  for (id in names(manifest$scenarios)) {
    entry <- sim_repair_scenario_entry(manifest$scenarios[[id]], id, paths)
    result_file <- entry$result_file
    checkpoint_file <- entry$checkpoint_file
    target <- entry$nsim_target

    if (sim_safe_file_exists(result_file)) {
      env <- new.env()
      load(result_file, envir = env)
      result <- env$result
      done <- sim_count_done_replicates(result$res.array)
      attempted <- sim_count_attempted_replicates(result$res.array, result$sim_fail)
      n_failed <- sum(!is.na(result$sim_fail))
      entry$nsim_done <- done
      entry$n_failed <- n_failed
      if (attempted >= target) {
        entry$status <- if (done > 0L) "completed" else "failed"
        entry$next_replicate <- if (done > 0L) target + 1L else 1L
        if (done > 0L && (is.null(entry$completed_at) || is.na(entry$completed_at))) {
          entry$completed_at <- as.character(file.mtime(result_file))
        }
        if (done > 0L) {
          sim_clear_checkpoint(id, paths)
        }
      } else {
        entry$status <- "partial"
        entry$next_replicate <- attempted + 1L
      }
    } else if (sim_safe_file_exists(checkpoint_file)) {
      ck <- readRDS(checkpoint_file)
      done <- ck$nsim_done %||% sim_count_done_replicates(ck$res.array)
      attempted <- sim_count_attempted_replicates(ck$res.array, ck$sim_fail)
      entry$nsim_done <- done
      entry$n_failed <- sum(!is.na(ck$sim_fail))
      entry$status <- if (attempted >= target) {
        if (done > 0L) "completed" else "failed"
      } else {
        "partial"
      }
      entry$next_replicate <- if (attempted >= target && done == 0L) {
        1L
      } else {
        ck$next_replicate %||% (attempted + 1L)
      }
      entry$started_at <- entry$started_at %||% ck$started_at
      entry$elapsed_sec <- ck$elapsed_sec
      entry$sec_per_replicate <- ck$sec_per_replicate
    } else {
      # No result or checkpoint on disk — reset unless currently mid-run elsewhere.
      if (entry$status != "running") {
        entry$status <- "pending"
        entry$nsim_done <- 0L
        entry$n_failed <- 0L
        entry$next_replicate <- 1L
        entry$completed_at <- NA_character_
      }
    }

    manifest$scenarios[[id]] <- entry
  }

  sim_save_manifest(manifest, paths)
}

sim_progress_update <- function(
    id,
    status = NULL,
    nsim_done = NULL,
    n_failed = NULL,
    next_replicate = NULL,
    started_at = NULL,
    completed_at = NULL,
    elapsed_sec = NULL,
    sec_per_replicate = NULL,
    manifest = NULL,
    paths = sim_paths()) {

  if (is.null(manifest)) {
    manifest <- sim_load_manifest(paths, create = TRUE)
  }
  if (!id %in% names(manifest$scenarios)) {
    stop("Unknown scenario id: ", id)
  }

  entry <- sim_repair_scenario_entry(manifest$scenarios[[id]], id, paths)
  if (!is.null(status)) entry$status <- status
  if (!is.null(nsim_done)) entry$nsim_done <- nsim_done
  if (!is.null(n_failed)) entry$n_failed <- n_failed
  if (!is.null(next_replicate)) entry$next_replicate <- next_replicate
  if (!is.null(started_at)) entry$started_at <- started_at
  if (!is.null(completed_at)) entry$completed_at <- completed_at
  if (!is.null(elapsed_sec)) entry$elapsed_sec <- elapsed_sec
  if (!is.null(sec_per_replicate)) entry$sec_per_replicate <- sec_per_replicate
  entry$updated_at <- sim_now_utc()
  manifest$scenarios[[id]] <- entry
  sim_save_manifest(manifest, paths)
  invisible(entry)
}

sim_progress_summary <- function(manifest = NULL, paths = sim_paths(), sync = TRUE) {
  if (sync) {
    manifest <- sim_sync_manifest_from_disk(manifest, paths)
  } else if (is.null(manifest)) {
    manifest <- sim_load_manifest(paths, create = TRUE)
  }

  entries <- manifest$scenarios
  statuses <- vapply(entries, function(x) x$status, character(1))
  done <- vapply(entries, function(x) x$nsim_done, integer(1))
  target <- vapply(entries, function(x) x$nsim_target, integer(1))
  failed <- vapply(entries, function(x) x$n_failed, integer(1))

  data.frame(
    scenario_id = names(entries),
    overlap = vapply(entries, function(x) x$overlap, character(1)),
    sample_size = vapply(entries, function(x) x$sample_size, character(1)),
    trial_type = vapply(entries, function(x) x$trial_type, character(1)),
    status = statuses,
    nsim_done = done,
    nsim_target = target,
    pct_done = round(100 * done / target, 1),
    n_failed = failed,
    next_replicate = vapply(entries, function(x) x$next_replicate, integer(1)),
    sec_per_rep = vapply(entries, function(x) x$sec_per_replicate %||% NA_real_, numeric(1)),
    stringsAsFactors = FALSE
  )
}

sim_format_duration <- function(seconds) {
  if (is.na(seconds) || !is.finite(seconds)) {
    return("—")
  }
  seconds <- as.integer(round(seconds))
  if (seconds < 60) {
    return(paste0(seconds, "s"))
  }
  if (seconds < 3600) {
    return(paste0(seconds %/% 60, "m ", seconds %% 60, "s"))
  }
  paste0(seconds %/% 3600, "h ", (seconds %% 3600) %/% 60, "m")
}

sim_print_progress <- function(paths = sim_paths(), sync = TRUE) {
  manifest <- if (sync) {
    sim_sync_manifest_from_disk(paths = paths)
  } else {
    sim_load_manifest(paths, create = TRUE)
  }

  tab <- sim_progress_summary(manifest, paths, sync = FALSE)
  n_completed <- sum(tab$status == "completed")
  n_partial <- sum(tab$status == "partial")
  n_running <- sum(tab$status == "running")
  n_pending <- sum(tab$status == "pending")
  n_total <- nrow(tab)

  cat("\nPSC simulation progress\n")
  cat(paste(rep("=", 72), collapse = ""), "\n", sep = "")
  cat("Manifest:", sim_manifest_path(paths), "\n")
  cat("Updated: ", manifest$updated_at, "\n\n", sep = "")
  cat(
    "Scenarios: ", n_completed, " completed | ",
    n_partial, " partial | ",
    n_running, " running | ",
    n_pending, " pending | ",
    n_total, " total\n",
    sep = ""
  )

  overall_done <- sum(tab$nsim_done)
  overall_target <- sum(tab$nsim_target)
  cat(
    "Replicates: ", overall_done, " / ", overall_target,
    " (", round(100 * overall_done / overall_target, 1), "%)\n\n",
    sep = ""
  )

  running <- tab[tab$status %in% c("running", "partial"), , drop = FALSE]
  if (nrow(running) > 0) {
    cat("Resume from:\n")
    for (i in seq_len(nrow(running))) {
      r <- running[i, ]
      eta <- NA_real_
      if (!is.na(r$sec_per_rep) && r$nsim_done < r$nsim_target) {
        eta <- r$sec_per_rep * (r$nsim_target - r$nsim_done)
      }
      cat(
        "  - ", r$scenario_id,
        "  replicate ", r$next_replicate, "/", r$nsim_target,
        if (!is.na(eta)) paste0("  (~", sim_format_duration(eta), " remaining)") else "",
        "\n",
        sep = ""
      )
    }
    cat("\n")
  }

  if (n_pending > 0) {
    next_pending <- tab[tab$status == "pending", , drop = FALSE][1, ]
    cat("Next pending scenario: ", next_pending$scenario_id, "\n\n", sep = "")
  }

  print(tab[, c("scenario_id", "status", "nsim_done", "nsim_target", "pct_done", "n_failed", "next_replicate")],
        row.names = FALSE)

  if (n_completed == n_total) {
    cat("\nAll scenarios complete.\n")
  } else {
    cat("\nTo resume: source('run_resume.R')  or  Rscript run_resume.R\n")
    cat("To check only: source('show_progress.R')\n")
  }

  invisible(tab)
}

sim_save_checkpoint <- function(
    id,
    res.array,
    sim_fail,
    nsim_done,
    next_replicate,
    nsim_target,
    started_at,
    elapsed_sec,
    sec_per_replicate,
    paths = sim_paths(),
    update_progress = TRUE) {

  dir.create(sim_checkpoint_dir(paths), recursive = TRUE, showWarnings = FALSE)
  ck <- list(
    scenario_id = id,
    res.array = res.array,
    sim_fail = sim_fail,
    nsim_done = nsim_done,
    next_replicate = next_replicate,
    nsim_target = nsim_target,
    started_at = started_at,
    elapsed_sec = elapsed_sec,
    sec_per_replicate = sec_per_replicate,
    saved_at = sim_now_utc()
  )
  saveRDS(ck, sim_checkpoint_path(id, paths))
  if (update_progress) {
    done_now <- as.integer(nsim_done %||% 0L)
    next_rep <- as.integer(next_replicate %||% 1L)
    target <- as.integer(nsim_target %||% length(sim_fail))
    status <- if (!is.na(next_rep) && !is.na(target) && next_rep > target) {
      if (done_now > 0L) "completed" else "failed"
    } else {
      "partial"
    }
    sim_progress_update(
      id,
      status = status,
      nsim_done = done_now,
      n_failed = sum(!is.na(sim_fail)),
      next_replicate = next_rep,
      started_at = started_at,
      elapsed_sec = elapsed_sec,
      sec_per_replicate = sec_per_replicate,
      paths = paths
    )
  }
  invisible(ck)
}

sim_load_checkpoint <- function(id, paths = sim_paths()) {
  ck_file <- sim_checkpoint_path(id, paths)
  if (!sim_safe_file_exists(ck_file)) {
    return(NULL)
  }
  readRDS(ck_file)
}

sim_clear_checkpoint <- function(id, paths = sim_paths()) {
  ck_file <- sim_checkpoint_path(id, paths)
  if (sim_safe_file_exists(ck_file)) {
    file.remove(ck_file)
  }
  invisible(TRUE)
}

sim_scenario_is_complete <- function(id, nsim = NULL, paths = sim_paths()) {
  state <- sim_load_resume_state(id, nsim = nsim, paths = paths)
  if (is.null(state)) {
    return(FALSE)
  }
  # All-failed runs must not count as complete (e.g. hybrid bug).
  state$attempted >= state$nsim_target && state$done > 0L
}

sim_combine_summaries <- function(paths = sim_paths()) {
  manifest <- sim_sync_manifest_from_disk(paths = paths)
  parts <- list()
  for (id in names(manifest$scenarios)) {
    if (!sim_scenario_is_complete(id, paths = paths)) {
      next
    }
    env <- new.env()
    load(sim_scenario_result_path(id, paths), envir = env)
    parts[[id]] <- cbind(scenario = id, env$result$summary)
  }
  if (length(parts) == 0) {
    return(NULL)
  }
  out <- do.call(rbind, parts)
  out_file <- file.path(paths$results_dir, "all_scenarios_summary.csv")
  write.csv(out, out_file, row.names = FALSE)
  cat("Combined summary (completed scenarios only):", out_file, "\n")
  invisible(out)
}
