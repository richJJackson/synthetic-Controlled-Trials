# 100 hybrid replicates on the two outcome-shifted populations.
# Moderate: dataset_2_shift, outcome_shift 0.20.
# Poor: dataset_3_shift, outcome_shift 0.45.
# Historical controls remain the 500 draws from dataset 1.
# Writes data/hybrid_shift_100.csv and does not touch hybrid_spline.csv.

Sys.setenv(HYBRID_SOURCE_ONLY = "1")
sim_dir <- normalizePath("~/Documents/GitHub/synthetic-Controlled-Trials/Examples/psc_sim", mustWork = TRUE)
source(file.path(sim_dir, "hybrid_100.R"))

truth <- log(0.7)
n_reps <- 200L
n_cores <- 4L
batch_reps <- 10L
overlaps <- c(
  dataset_2_shift = "Moderate overlap, outcome shift",
  dataset_3_shift = "Poor overlap, outcome shift"
)
seeds <- c(dataset_2_shift = 14L, dataset_3_shift = 15L)
out_path <- file.path(sim_dir, "data/hybrid_shift_100.csv")
done_df <- if (file.exists(out_path)) utils::read.csv(out_path) else NULL

run_one <- function(i) {
  rep_i <- tasks$rep_i[i]
  nm <- tasks$nm[i]
  seed_h <- unname(seeds[[nm]]) + (rep_i - 1L) * 1000L
  one <- NULL
  utils::capture.output({
    one <- run_rep(seed_h, nm)
  })
  one$rep <- rep_i
  one$overlap <- unname(overlaps[[nm]])
  one
}

cover <- function(est, lo, hi) {
  use <- is.finite(est)
  c(
    n = sum(use),
    mean = mean(est[use]),
    median = stats::median(est[use]),
    bias = mean(est[use] - truth),
    cover = sum(lo[use] <= truth & hi[use] >= truth),
    acil = mean((hi - lo)[use])
  )
}

report <- function(tab) {
  for (ov in overlaps) {
    sub <- tab[tab$overlap == ov, ]
    summary_tbl <- rbind(
      pooled = cover(sub$est_pooled, sub$lo_pooled, sub$hi_pooled),
      population_sc = cover(sub$est_sc, sub$lo_sc, sub$hi_sc),
      case_weighted = cover(sub$est_cw, sub$lo_cw, sub$hi_cw),
      psc = cover(sub$est_psc, sub$lo_psc, sub$hi_psc),
      commensurate = cover(sub$est_comm, sub$lo_comm, sub$hi_comm),
      leap = cover(sub$est_leap, sub$lo_leap, sub$hi_leap)
    )
    cat(ov, "replicates", nrow(sub), "PSC accepted", sum(sub$psc_accepted), "\n")
    print(round(summary_tbl, 3))
    psc_use <- is.finite(sub$est_psc)
    cat(
      "tau", round(mean(sub$tau), 2),
      "sum_w", round(mean(sub$sum_w[psc_use]), 1),
      "n_exchangeable", round(mean(sub$n_exchangeable), 1),
      "\n"
    )
  }
}

next_rep <- if (is.null(done_df)) 1L else max(done_df$rep) + 1L
while (next_rep <= n_reps) {
  rep_to <- min(next_rep + batch_reps - 1L, n_reps)
  tasks <- expand.grid(
    rep_i = next_rep:rep_to,
    nm = names(overlaps),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  t0 <- proc.time()[["elapsed"]]
  new_rows <- parallel::mclapply(seq_len(nrow(tasks)), run_one, mc.cores = n_cores)
  new_df <- do.call(rbind, new_rows)
  done_df <- if (is.null(done_df)) new_df else rbind(done_df, new_df[, names(done_df), drop = FALSE])
  utils::write.csv(done_df, out_path, row.names = FALSE)
  cat("UPDATE through replicate", rep_to, "sec", round(proc.time()[["elapsed"]] - t0, 1), "\n")
  report(done_df)
  flush.console()
  next_rep <- rep_to + 1L
}
cat("Truth", truth, "\n")
cat("Saved", out_path, "\n")
