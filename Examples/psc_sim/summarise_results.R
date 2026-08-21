# Summarise completed PSC simulation scenarios.
#
# Usage (from psc_sim/):
#   source("summarise_results.R")
#
# Or one scenario only:
#   source("summarise_results.R"); summarise_scenario("no_overlap_small_single_arm")
#
# Or from the shell:
#   Rscript summarise_results.R

if (!exists("PSC_SIM_DIR")) {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) {
    PSC_SIM_DIR <- normalizePath(dirname(sub("^--file=", "", file_arg[1])), mustWork = TRUE)
  } else {
    PSC_SIM_DIR <- normalizePath(
      "~/Documents/GitHub/synthetic-Controlled-Trials/Examples/psc_sim"
    )
  }
}
assign("PSC_SIM_DIR", PSC_SIM_DIR, envir = .GlobalEnv)
Sys.setenv(PSC_SIM_DIR = PSC_SIM_DIR)

source(file.path(PSC_SIM_DIR, "setup.R"))
source(file.path(PSC_SIM_DIR, "config.R"))
source(file.path(PSC_SIM_DIR, "summarise.R"))
source(file.path(PSC_SIM_DIR, "progress.R"))

sim_round_summary <- function(x, digits = 4) {
  num <- vapply(x, is.numeric, logical(1))
  x[num] <- lapply(x[num], function(col) round(col, digits))
  x
}

print_sim_summary <- function(summary_tab, result, title = NULL) {
  title <- title %||% result$scenario_id
  cat("\n", title, "\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
  cat(
    "Overlap: ", result$overlap,
    "  |  Sample size: ", result$sample_size,
    "  |  Trial type: ", result$trial_type, "\n",
    "n_cont = ", result$n_cont, ", n_trt = ", result$n_trt,
    ", nsim = ", result$nsim,
    ", failed = ", sum(!is.na(result$sim_fail)), "\n",
    "True log(HR) = ", round(result$true_beta, 4),
    "  (HR = ", round(exp(result$true_beta), 4), ")\n",
    "Coverage uses stored interval limits (95% HPD for PSC; Wald for others).\n\n",
    sep = ""
  )
  print(sim_round_summary(summary_tab), row.names = FALSE)
  invisible(summary_tab)
}

load_scenario_result <- function(scenario_id, paths = sim_paths()) {
  result_file <- sim_scenario_result_path(scenario_id, paths)
  if (!file.exists(result_file)) {
    stop("Result file not found: ", result_file)
  }
  env <- new.env()
  load(result_file, envir = env)
  env$result
}

completed_scenario_ids <- function(paths = sim_paths()) {
  manifest <- sim_sync_manifest_from_disk(paths = paths)
  ids <- names(manifest$scenarios)[
    vapply(manifest$scenarios, function(x) x$status == "completed", logical(1))
  ]
  if (length(ids) == 0) {
    result_files <- list.files(
      paths$results_dir,
      pattern = "\\.RData$",
      full.names = FALSE
    )
    ids <- sub("\\.RData$", "", result_files)
  }
  sort(ids)
}

summarise_scenario <- function(
    scenario_id,
    paths = sim_paths(),
    save_csv = TRUE,
    plot = TRUE) {

  result <- load_scenario_result(scenario_id, paths)
  summary_tab <- result$summary
  if (is.null(summary_tab)) {
    method_names <- if (result$trial_type == "hybrid") {
      HYBRID_METHOD_NAMES
    } else {
      SA_METHOD_NAMES
    }
    summary_tab <- summarise_res_array(
      result$res.array,
      true_val = result$true_beta,
      method_names = method_names
    )
  }

  print_sim_summary(summary_tab, result)

  if (save_csv) {
    out_file <- file.path(
      paths$results_dir,
      paste0(scenario_id, "_summary.csv")
    )
    write.csv(summary_tab, out_file, row.names = FALSE)
    cat("\nSaved:", out_file, "\n")
  }

  if (plot) {
    plot_file <- file.path(
      paths$results_dir,
      paste0(scenario_id, "_bias.png")
    )
    plot_scenario_summary(summary_tab, result, plot_file)
    cat("Plot saved:", plot_file, "\n")
  }

  invisible(list(result = result, summary = summary_tab))
}

plot_scenario_summary <- function(summary_tab, result, outfile = NULL) {
  if (!is.null(outfile)) {
    grDevices::png(outfile, width = 1000, height = 700, res = 120)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mar = c(10, 4, 3, 1))
  bias <- summary_tab$mean_bias
  cols <- ifelse(abs(bias) < 0.05, "steelblue", "coral3")
  barplot(
    bias,
    names.arg = summary_tab$method,
    las = 2,
    col = cols,
    ylab = expression("Mean bias (" * log * " HR)"),
    main = paste0("Bias by method: ", result$scenario_id),
    cex.names = 0.75
  )
  abline(h = 0, lty = 2)
}

summarise_all_completed <- function(
    paths = sim_paths(),
    save_csv = TRUE,
    plot = FALSE) {

  ids <- completed_scenario_ids(paths)
  if (!length(ids)) {
    message("No completed scenarios found in ", paths$results_dir)
    return(invisible(NULL))
  }

  cat("Completed scenarios (", length(ids), "):\n", sep = "")
  cat("  ", paste(ids, collapse = "\n   "), "\n\n", sep = "")

  parts <- lapply(ids, function(id) {
    out <- summarise_scenario(id, paths = paths, save_csv = FALSE, plot = plot)
    cbind(scenario_id = id, out$summary)
  })

  combined <- do.call(rbind, parts)
  rownames(combined) <- NULL

  if (save_csv) {
    per_scenario <- file.path(paths$results_dir, paste0(ids, "_summary.csv"))
    for (i in seq_along(ids)) {
      write.csv(parts[[i]][, setdiff(names(parts[[i]]), "scenario_id")],
                per_scenario[i], row.names = FALSE)
    }
    combined_file <- file.path(paths$results_dir, "all_scenarios_summary.csv")
    write.csv(combined, combined_file, row.names = FALSE)
    cat("\nCombined summary:", combined_file, "\n")
  }

  invisible(combined)
}

CORE_METHOD_NAMES <- c(
  "Pooled (unadj)",
  "PSC full",
  "SC full",
  "Bayes hist (info)",
  "Bayes hist (case-weighted)"
)

CORE_METHOD_LABELS <- c(
  "Pooled (unadj)" = "Pooled",
  "PSC full" = "PSC",
  "SC full" = "SC",
  "Bayes hist (info)" = "Bayes (info)",
  "Bayes hist (case-weighted)" = "Bayes (CW)"
)

prepare_core_summary <- function(paths = sim_paths()) {
  combined_file <- file.path(paths$results_dir, "all_scenarios_summary.csv")
  if (!file.exists(combined_file)) {
    combined <- summarise_all_completed(paths = paths, save_csv = TRUE, plot = FALSE)
  } else {
    combined <- read.csv(combined_file, stringsAsFactors = FALSE)
  }
  combined <- combined[combined$method %in% CORE_METHOD_NAMES, , drop = FALSE]
  if (!nrow(combined)) {
    stop("No rows for core methods in summary data.")
  }

  id_col <- if ("scenario_id" %in% names(combined)) "scenario_id" else "scenario"
  parsed <- lapply(combined[[id_col]], parse_scenario_id)
  combined$overlap <- vapply(parsed, `[[`, "", "overlap")
  combined$sample_size <- vapply(parsed, `[[`, "", "sample_size")
  combined$trial_type <- vapply(parsed, `[[`, "", "trial_type")
  combined <- combined[combined$trial_type == "single_arm", , drop = FALSE]
  combined$method_label <- CORE_METHOD_LABELS[combined$method]

  overlap_lev <- intersect(
    c("no_overlap", "moderate_overlap", "large_overlap"),
    unique(combined$overlap)
  )
  size_lev <- intersect(c("small", "medium", "large"), unique(combined$sample_size))
  combined$overlap <- factor(combined$overlap, levels = overlap_lev)
  combined$sample_size <- factor(
    combined$sample_size,
    levels = size_lev,
    labels = paste0(toupper(substring(size_lev, 1, 1)), substring(size_lev, 2))
  )
  combined$method_label <- factor(
    combined$method_label,
    levels = unname(CORE_METHOD_LABELS[CORE_METHOD_NAMES])
  )
  combined
}

plot_core_methods_overview <- function(
    combined = NULL,
    paths = sim_paths(),
    outfile = file.path(paths$results_dir, "core_methods_overview.png"),
    width = 11,
    height = 7,
    dpi = 150) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_core_methods_overview().")
  }

  if (is.null(combined)) {
    combined <- prepare_core_summary(paths)
  }

  overlap_labels <- c(
    no_overlap = "No overlap",
    moderate_overlap = "Moderate overlap",
    large_overlap = "Large overlap"
  )
  combined$overlap_label <- overlap_labels[as.character(combined$overlap)]

  pal <- c(
    "Pooled" = "#4d4d4d",
    "PSC" = "#1b7837",
    "SC" = "#762a83",
    "Bayes (info)" = "#d73027",
    "Bayes (CW)" = "#e08214"
  )

  p <- ggplot2::ggplot(
    combined,
    ggplot2::aes(x = mean_bias, y = coverage, colour = method_label)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey70") +
    ggplot2::geom_point(size = 2.8, alpha = 0.9) +
    ggplot2::facet_grid(
      sample_size ~ overlap_label,
      labeller = ggplot2::label_value
    ) +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::scale_x_continuous(
      name = expression("Bias (" * log * " HR)"),
      limits = c(-0.9, 0.15)
    ) +
    ggplot2::scale_y_continuous(
      name = "Coverage of 95% intervals",
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      title = "Simulation results: single-arm scenarios",
      subtitle = paste(
        "True log(HR) =", round(SIM_DGP$beta, 3),
        "(HR =", round(exp(SIM_DGP$beta), 2), ")"
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box.spacing = ggplot2::unit(0.2, "cm"),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(outfile, p, width = width, height = height, dpi = dpi, bg = "white")
  invisible(list(plot = p, file = outfile, data = combined))
}

plot_core_methods_acil <- function(
    combined = NULL,
    paths = sim_paths(),
    outfile = file.path(paths$results_dir, "core_methods_acil.png"),
    width = 11,
    height = 6.5,
    dpi = 150) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_core_methods_acil().")
  }

  if (is.null(combined)) {
    combined <- prepare_core_summary(paths)
  }

  overlap_labels <- c(
    no_overlap = "No overlap",
    moderate_overlap = "Moderate overlap",
    large_overlap = "Large overlap"
  )
  combined$overlap_label <- overlap_labels[as.character(combined$overlap)]

  pal <- c(
    "Pooled" = "#4d4d4d",
    "PSC" = "#1b7837",
    "SC" = "#762a83",
    "Bayes (info)" = "#d73027",
    "Bayes (CW)" = "#e08214"
  )

  p <- ggplot2::ggplot(
    combined,
    ggplot2::aes(x = method_label, y = acil, fill = method_label)
  ) +
    ggplot2::geom_col(width = 0.72, show.legend = FALSE) +
    ggplot2::facet_grid(
      sample_size ~ overlap_label,
      labeller = ggplot2::label_value
    ) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::scale_y_continuous(
      name = "Average confidence interval length (log HR)",
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      x = NULL,
      title = "Average confidence interval length",
      subtitle = "Same scenarios and methods as bias/coverage figure"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(outfile, p, width = width, height = height, dpi = dpi, bg = "white")
  invisible(list(plot = p, file = outfile, data = combined))
}

# --- run when sourced or executed as a script ---
ids <- completed_scenario_ids()
if (length(ids) == 0) {
  message("No completed scenario results yet. Check progress with source('show_progress.R').")
} else if (length(ids) == 1) {
  summarise_scenario(ids[[1]])
} else {
  summarise_all_completed()
}
