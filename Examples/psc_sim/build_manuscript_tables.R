# Build manuscript tables and coverage-vs-bias figures for core methods.
#
# Usage (from psc_sim/):
#   Rscript build_manuscript_tables.R

Sys.setenv(PSC_SIM_DIR = normalizePath("."))
source("config.R")
suppressPackageStartupMessages(library(ggplot2))

paths <- sim_paths()
d <- read.csv(file.path(paths$results_dir, "all_scenarios_summary.csv"), stringsAsFactors = FALSE)

minus <- "\u2212"  # unicode minus

fmt_num <- function(x, digits = 3) {
  sprintf(paste0("%s%0.", digits, "f"), ifelse(x < 0, minus, ""), abs(x))
}

fmt_mean_se <- function(est, se) {
  paste0(fmt_num(est, 3), " (", sprintf("%.3f", se), ")")
}

fmt_cov <- function(x) sprintf("%.1f", 100 * x)

overlap_lab <- c(
  no_overlap = "Small",
  moderate_overlap = "Medium",
  large_overlap = "Large"
)
size_lab <- c(small = "Small", medium = "Medium", large = "Large")

parse_row <- function(id) {
  p <- parse_scenario_id(id)
  data.frame(
    overlap_code = p$overlap,
    size_code = p$sample_size,
    Overlap = unname(overlap_lab[[p$overlap]]),
    `Sample size` = unname(size_lab[[p$sample_size]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

build_block <- function(d, trial_suffix, method_map, design_label) {
  ids <- d$scenario_id[grepl(paste0("_", trial_suffix, "$"), d$scenario_id)]
  ids <- unique(ids)
  parsed <- do.call(rbind, lapply(ids, parse_row))
  parsed$scenario_id <- ids
  out <- list()
  size_ord <- c("Small", "Medium", "Large")
  ov_ord <- c("Small", "Medium", "Large")
  for (disp in names(method_map)) {
    raw <- method_map[[disp]]
    for (sz in size_ord) {
      for (ov in ov_ord) {
        sid <- parsed$scenario_id[parsed$`Sample size` == sz & parsed$Overlap == ov]
        r <- d[d$scenario_id == sid & d$method == raw, , drop = FALSE]
        if (!nrow(r)) {
          stop("Missing: ", raw, " in ", sid)
        }
        out[[length(out) + 1L]] <- data.frame(
          Design = design_label,
          Method = disp,
          `Sample size` = sz,
          Overlap = ov,
          `Mean (se)` = fmt_mean_se(r$mean_est, r$mean_se),
          ACIL = fmt_num(r$acil, 3),
          Coverage = fmt_cov(r$coverage),
          Bias = fmt_num(r$mean_bias, 3),
          scenario_id = sid,
          method_raw = raw,
          mean_est = r$mean_est,
          mean_se = r$mean_se,
          acil = r$acil,
          coverage = r$coverage,
          mean_bias = r$mean_bias,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    }
  }
  do.call(rbind, out)
}

sa_map <- list(
  "Unadjusted Cox" = "Pooled (unadj)",
  "Personalised Synthetic Controls" = "PSC full",
  "Synthetic Controls" = "SC full",
  "Bayesian Case-Weighted" = "Bayes hist (case-weighted)"
)
hy_map <- list(
  "Unadjusted Cox" = "RCT (unadj)",
  "Personalised Synthetic Controls" = "PSC full (combined)",
  "Synthetic Controls" = "SC full",
  "Bayesian Case-Weighted" = "Bayes case-weighted",
  "Bayesian Commensurate Priors" = "Bayes commensurate"
)

sa <- build_block(d, "single_arm", sa_map, "Fully synthetic")
hy <- build_block(d, "hybrid", hy_map, "Partially synthetic")

sa_out <- sa[, c("Method", "Sample size", "Overlap", "Mean (se)", "ACIL", "Coverage", "Bias")]
hy_out <- hy[, c("Method", "Sample size", "Overlap", "Mean (se)", "ACIL", "Coverage", "Bias")]

man_dir <- file.path(dirname(dirname(paths$sim_dir)), "Manuscript")
dir.create(man_dir, showWarnings = FALSE)
sa_csv <- file.path(man_dir, "Table_fully_synthetic_core.csv")
hy_csv <- file.path(man_dir, "Table_hybrid_core.csv")
both_csv <- file.path(man_dir, "Table_core_methods_all.csv")
write.csv(sa_out, sa_csv, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(hy_out, hy_csv, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(
  rbind(sa, hy)[, c("Design", "Method", "Sample size", "Overlap", "Mean (se)", "ACIL", "Coverage", "Bias")],
  both_csv, row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(sa_out, file.path(paths$results_dir, "fully_synthetic_core_table.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(hy_out, file.path(paths$results_dir, "hybrid_core_table.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("Saved tables:\n ", sa_csv, "\n ", hy_csv, "\n")

# ---- coverage vs bias figures ----
# Overlap of nearby scenario summaries hid points: filled pch 15–17 with
# alpha encoding, a 0–100% y-axis, and several same-method/same-size pairs
# within ~0.007 bias / ~0.001 coverage (hybrid PSC and commensurate priors).
# Use outlined shapes; small overlap is hollow, large is solid; autoscale y.
pt_shapes <- c(Small = 21, Medium = 24, Large = 22)

prepare_cov_bias <- function(dat, method_levels = NULL) {
  dat$Overlap <- factor(dat$Overlap, levels = c("Small", "Medium", "Large"))
  dat$`Sample size` <- factor(dat$`Sample size`, levels = c("Small", "Medium", "Large"))
  if (is.null(method_levels)) {
    method_levels <- unique(dat$Method)
  }
  dat$Method <- factor(dat$Method, levels = method_levels)
  # Deterministic offsets (not random jitter): sample size along bias,
  # overlap along coverage. Separates stacked same-method points such as
  # hybrid commensurate-prior medium size (bias 0.023 vs 0.030, coverage
  # 0.922 vs 0.921) and fully synthetic case-weighted large overlap
  # (three sizes at ~0.05 bias, ~97.6% coverage).
  sz <- c(Small = -1, Medium = 0, Large = 1)[as.character(dat$`Sample size`)]
  ov <- c(Small = 1, Medium = 0, Large = -1)[as.character(dat$Overlap)]
  dat$bias_plot <- dat$mean_bias + 0.0048 * sz
  dat$cov_plot <- dat$coverage + 0.0030 * ov
  dat[order(dat$Overlap, decreasing = TRUE), , drop = FALSE]
}

cov_bias_plot <- function(dat, title, subtitle, method_cols, facet = FALSE) {
  dat <- prepare_cov_bias(dat, method_levels = names(method_cols))
  p <- ggplot(
    dat,
    aes(
      x = bias_plot, y = cov_plot,
      colour = Method, fill = Method, shape = `Sample size`
    )
  ) +
    geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50", linewidth = 0.4) +
    geom_point(aes(alpha = Overlap), size = 3.0, stroke = 0) +
    geom_point(alpha = 1, size = 3.0, stroke = 0.7, fill = NA) +
    scale_colour_manual(values = method_cols, name = "Method") +
    scale_fill_manual(values = method_cols, name = "Method") +
    scale_shape_manual(values = pt_shapes, name = "Sample size") +
    scale_alpha_manual(
      values = c(Small = 0, Medium = 0.55, Large = 1),
      name = "Overlap"
    ) +
    scale_y_continuous(
      labels = function(x) paste0(format(100 * x, nsmall = 0), "%"),
      expand = expansion(mult = c(0.08, 0.06))
    ) +
    scale_x_continuous(expand = expansion(mult = 0.06)) +
    labs(
      x = "Mean bias (log hazard ratio)",
      y = "Coverage of nominal 95% interval",
      title = title,
      subtitle = subtitle
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, colour = "grey30")
    ) +
    guides(
      colour = guide_legend(
        order = 1, nrow = 2,
        override.aes = list(
          shape = 21, alpha = 1, stroke = 0.7, size = 3,
          fill = unname(method_cols)
        )
      ),
      fill = "none",
      shape = guide_legend(
        order = 2, nrow = 1,
        override.aes = list(fill = "grey55", colour = "grey20", alpha = 1, stroke = 0.7)
      ),
      alpha = guide_legend(
        order = 3, nrow = 1,
        override.aes = list(
          shape = 21, size = 3, fill = "grey35", colour = "grey20", stroke = 0.7
        )
      )
    )
  if (facet) {
    p <- p + facet_wrap(~Design, nrow = 1, scales = "free_y")
  }
  p
}

plot_fun <- function(dat, title, subtitle, outfile_png, outfile_pdf, method_cols) {
  p <- cov_bias_plot(dat, title, subtitle, method_cols, facet = FALSE)
  ggsave(outfile_png, p, width = 8.2, height = 6.2, dpi = 300)
  ggsave(outfile_pdf, p, width = 8.2, height = 6.2)
  invisible(p)
}

sub <- paste0(
  "Nine points per method (3 sample sizes × 3 overlap levels); true log HR = ",
  round(SIM_DGP$beta, 3), " (HR = ", round(exp(SIM_DGP$beta), 2), "). ",
  "Fill: open = small overlap, solid = large. Points are offset slightly ",
  "by sample size (horizontal) and overlap (vertical) to reduce overplotting."
)

# Figures omit Unadjusted Cox; tables above retain it.
sa_fig <- sa[sa$Method != "Unadjusted Cox", , drop = FALSE]
hy_fig <- hy[hy$Method != "Unadjusted Cox", , drop = FALSE]

sa_cols <- c(
  "Personalised Synthetic Controls" = "#D55E00",
  "Synthetic Controls" = "#CC79A7",
  "Bayesian Case-Weighted" = "#009E73"
)
hy_cols <- c(
  "Personalised Synthetic Controls" = "#D55E00",
  "Synthetic Controls" = "#CC79A7",
  "Bayesian Case-Weighted" = "#009E73",
  "Bayesian Commensurate Priors" = "#E69F00"
)

sa_png <- file.path(paths$results_dir, "coverage_vs_bias_fully_synthetic.png")
sa_pdf <- file.path(paths$results_dir, "coverage_vs_bias_fully_synthetic.pdf")
hy_png <- file.path(paths$results_dir, "coverage_vs_bias_hybrid.png")
hy_pdf <- file.path(paths$results_dir, "coverage_vs_bias_hybrid.pdf")

p_sa <- plot_fun(sa_fig, "Fully synthetic trials", sub, sa_png, sa_pdf, sa_cols)
p_hy <- plot_fun(hy_fig, "Partially synthetic (hybrid) trials", sub, hy_png, hy_pdf, hy_cols)

# Combined two-panel Figure 1
both <- rbind(sa_fig, hy_fig)
both$Design <- factor(
  both$Design,
  levels = c("Fully synthetic", "Partially synthetic")
)
comb_cols <- hy_cols
comb_sub <- paste0(
  "Nine points per method (3 sample sizes × 3 overlap levels); true log HR = ",
  round(SIM_DGP$beta, 3), ".\n",
  "Open fill = small overlap, solid = large. Separate coverage scales; ",
  "slight offsets reduce overplotting. Commensurate priors: hybrid only."
)
p_comb <- cov_bias_plot(
  both,
  "Coverage versus bias across simulation scenarios",
  comb_sub,
  comb_cols,
  facet = TRUE
)

comb_png <- file.path(paths$results_dir, "coverage_vs_bias_sa_ra.png")
comb_pdf <- file.path(paths$results_dir, "coverage_vs_bias_sa_ra.pdf")
ggsave(comb_png, p_comb, width = 10.5, height = 6.2, dpi = 300)
ggsave(comb_pdf, p_comb, width = 10.5, height = 6.2)
write.csv(
  both[, c("Design", "Method", "Sample size", "Overlap", "mean_est", "mean_se",
           "acil", "coverage", "mean_bias", "scenario_id", "method_raw")],
  file.path(paths$results_dir, "coverage_vs_bias_sa_ra_data.csv"),
  row.names = FALSE
)

cat("Saved figures:\n ", sa_png, "\n ", hy_png, "\n ", comb_png, "\n")
cat("Figure methods:", paste(unique(as.character(both$Method)), collapse = "; "), "\n")
cat("SA figure rows:", nrow(sa_fig), " hybrid figure rows:", nrow(hy_fig), "\n")
