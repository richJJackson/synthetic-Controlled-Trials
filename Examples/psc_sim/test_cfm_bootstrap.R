# Quick test: bootstrap CFM coefficient draws vs MVN (vcov) in PSC.
#
# Usage (from psc_sim/):
#   Rscript test_cfm_bootstrap.R

source("config.R")
source("setup.R")
sim_source_all()
sim_load_packages()
suppressPackageStartupMessages(library(posterior))

paths <- sim_paths()
load(file.path(paths$data_dir, "populations.RData"))

overlap <- "moderate_overlap"
sample_size <- "medium"
seed <- 21320L
n_boot <- 200L
nsim <- 5000L
burn <- 500L
thin <- 2L

sizes <- SIM_SAMPLE_SIZES[[sample_size]]
cohorts <- sim_sample_cohorts(
  populations, overlap, sizes$n_cont, sizes$n_trt, seed = seed
)
prep <- sim_prepare_single_arm(cohorts$contDat, cohorts$trtDat)
aligned <- sim_factor_align(prep$contDat, prep$trtDat)
contDat <- aligned$contDat
trtDat <- aligned$trtDat

cat("Building CFM + bootstrap bank (B =", n_boot, ")...\n")
t_boot <- proc.time()[["elapsed"]]
cfm <- flexsurv::flexsurvspline(
  s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5,
  data = contDat,
  k = 3
)
cfmw <- as_pscCFM_dev(cfm)
co_names <- names(cfmw$co)
sig_mvn <- as.matrix(cfmw$sig)

boot_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(co_names))
colnames(boot_mat) <- co_names
n_ok <- 0L
attempts <- 0L
set.seed(seed + 99L)
while (n_ok < n_boot && attempts < n_boot * 5L) {
  attempts <- attempts + 1L
  idx <- sample.int(nrow(contDat), replace = TRUE)
  fit_b <- tryCatch(
    flexsurv::flexsurvspline(
      s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5,
      data = contDat[idx, , drop = FALSE],
      k = 3
    ),
    error = function(e) NULL
  )
  if (is.null(fit_b)) {
    next
  }
  b <- fit_b$coefficients
  if (!all(co_names %in% names(b))) {
    next
  }
  n_ok <- n_ok + 1L
  boot_mat[n_ok, ] <- as.numeric(b[co_names])
}
boot_mat <- boot_mat[seq_len(n_ok), , drop = FALSE]
cat(
  "Bootstrap fits OK:", n_ok, "/", attempts,
  " in ", round(proc.time()[["elapsed"]] - t_boot, 1), "s\n",
  sep = ""
)

boot_cov <- stats::cov(boot_mat)
mvn_var <- diag(sig_mvn)
boot_var <- diag(boot_cov)
cat("\nCFM coefficient SDs (MLE vcov vs bootstrap):\n")
print(round(cbind(mvn_sd = sqrt(mvn_var), boot_sd = sqrt(boot_var),
                  ratio = sqrt(boot_var) / sqrt(mvn_var)), 3))

#' Run pscEst with a chosen CFM draw mechanism.
run_psc_cfm <- function(cfm_draws = c("mvn", "bootstrap", "fixed"),
                        boot_bank = NULL) {
  cfm_draws <- match.arg(cfm_draws)
  pscOb <- pscData(cfmw, trtDat)
  pscOb <- init(pscOb)
  pscOb <- pscEst_start(
    pscOb,
    nsim = nsim,
    nchain = 1L,
    fixed_cfm = (cfm_draws == "fixed")
  )
  if (cfm_draws == "bootstrap") {
    stopifnot(!is.null(boot_bank), nrow(boot_bank) >= 20L)
    bank <- boot_bank
    pscOb$cfmPost <- function(x) {
      i <- sample.int(nrow(bank), size = x, replace = TRUE)
      if (x == 1L) {
        as.numeric(bank[i, ])
      } else {
        bank[i, , drop = FALSE]
      }
    }
    pscOb$fixed_cfm <- FALSE
  } else if (cfm_draws == "mvn") {
    # default cfmPost already MVN(co, sig)
    pscOb$fixed_cfm <- FALSE
  }
  pscOb <- pscEst_run(pscOb, nsim = nsim, nchain = 1L)
  pscOb <- postSummary(pscOb, thin = thin, burn = burn)
  beta <- as.numeric(as_draws(pscOb$draws)$beta_1)
  list(
    method = cfm_draws,
    beta = beta,
    median = stats::median(beta),
    sd = stats::sd(beta),
    hpd = sim_hpd_interval(beta),
    ess = as.numeric(posterior::ess_bulk(posterior::as_draws_matrix(
      data.frame(beta_1 = beta)
   ))),
    n = length(beta)
  )
}

cat("\nRunning PSC: fixed / MVN / bootstrap CFM (nsim =", nsim, ")...\n")
res_fixed <- run_psc_cfm("fixed")
res_mvn <- run_psc_cfm("mvn")
res_boot <- run_psc_cfm("bootstrap", boot_bank = boot_mat)

summ <- rbind(
  data.frame(
    cfm_draws = c("fixed", "mvn", "bootstrap"),
    n = c(res_fixed$n, res_mvn$n, res_boot$n),
    median = c(res_fixed$median, res_mvn$median, res_boot$median),
    sd = c(res_fixed$sd, res_mvn$sd, res_boot$sd),
    hpd_lo = c(res_fixed$hpd[["lower"]], res_mvn$hpd[["lower"]], res_boot$hpd[["lower"]]),
    hpd_hi = c(res_fixed$hpd[["upper"]], res_mvn$hpd[["upper"]], res_boot$hpd[["upper"]]),
    width = c(
      diff(res_fixed$hpd), diff(res_mvn$hpd), diff(res_boot$hpd)
    ),
    ess = c(res_fixed$ess, res_mvn$ess, res_boot$ess),
    stringsAsFactors = FALSE
  )
)
summ$covers_truth <- SIM_DGP$beta >= summ$hpd_lo & SIM_DGP$beta <= summ$hpd_hi

out_dir <- file.path(paths$results_dir, "gated")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(out_dir, "cfm_bootstrap_example.csv")
write.csv(summ, out_csv, row.names = FALSE)

cat("\nTrue log HR =", round(SIM_DGP$beta, 4), "\n")
cat("\nPSC beta_1 posterior summary:\n")
print(summ, row.names = FALSE)
cat("\nSD inflation bootstrap vs MVN:", round(res_boot$sd / res_mvn$sd, 3), "\n")
cat("Width inflation bootstrap vs MVN:", round(diff(res_boot$hpd) / diff(res_mvn$hpd), 3), "\n")
cat("\nSaved:", out_csv, "\n")
