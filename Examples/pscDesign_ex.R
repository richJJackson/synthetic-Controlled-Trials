### Comparing, SC, PSC and pooled analysis

###

install.packages("rstan")
library(ebal)
library(survival)
library(psc)
library(flexsurv)
library(rstan)

setwd("/Users/richardjackson/Documents/GitHub/psc")
devtools::load_all()


######## CONTROL POPULATION ############
### Creating population of control and experimental patients

N <- 5000

gam_0 <- 0.05
gam_1 <- log(0.6)
gam_2 <- log(0.8)
gam_3 <- log(0.6)
gam_4 <- log(1.2)
gam_5 <- log(0.9)
beta <- log(0.7)

X_1 <- rnorm(N)
X_2 <- rnorm(N)
X_3 <- rbinom(N, 1, 0.3)
X_4 <- rbinom(N, 1, 0.4)
X_5 <- rbinom(N, 1, 0.5)

lp <- exp(
  gam_0 +
    gam_1 * X_1 +
    gam_2 * X_2 +
    gam_3 * X_3 +
    gam_4 * X_4 +
    gam_5 * X_5 +
    beta
)
lpc <- exp(
  gam_0 +
    gam_1 * X_1 +
    gam_2 * X_2 +
    gam_3 * X_3 +
    gam_4 * X_4 +
    gam_5 * X_5
)

tm <- rexp(N, lp)
tmc <- rexp(N, lpc)
centm <- rexp(N, rep(0.02, N))
centm

### Creating dataset
cont <- data.frame(X_1, X_2, X_3, X_4, X_5, tm, tmc, centm)


######## Treated POPULATION ############
### Creating population of control and experimental patients
N <- 5000

gam_0 <- 0.05
gam_1 <- log(0.6)
gam_2 <- log(0.8)
gam_3 <- log(0.6)
gam_4 <- log(1.2)
gam_5 <- log(0.9)
beta <- log(0.7)

X_1 <- rnorm(N,1.5,1)
X_2 <- rnorm(N,0.5,1)
X_3 <- rbinom(N, 1, 0.6)
X_4 <- rbinom(N, 1, 0.7)
X_5 <- rbinom(N, 1, 0.8)


lp <- exp(
  gam_0 +
    gam_1 * X_1 +
    gam_2 * X_2 +
    gam_3 * X_3 +
    gam_4 * X_4 +
    gam_5 * X_5 +
    beta
)
lpc <- exp(
  gam_0 +
    gam_1 * X_1 +
    gam_2 * X_2 +
    gam_3 * X_3 +
    gam_4 * X_4 +
    gam_5 * X_5
)

tm <- rexp(N, lp)
tmc <- rexp(N, lpc)
centm <- rexp(N, rep(0.02, N))

### Creating dataset
treated <- data.frame(X_1, X_2, X_3, X_4, X_5, tm, tmc, centm)


#### Consoring patterns
cont$cen <- as.numeric(!(cont$centm < cont$tm))
cont$tm[which(cont$cen == 0)] <- cont$centm[which(cont$cen == 0)]
cont$cenc <- as.numeric(!(cont$centm < cont$tmc))
cont$tmc[which(cont$cenc == 0)] <- cont$centm[which(cont$cenc == 0)]

treated$cen <- as.numeric(!(treated$centm < treated$tm))
treated$cen
treated$tm[which(treated$cen == 0)] <- treated$centm[which(treated$cen == 0)]
treated$cen <- as.numeric(!(6 < treated$tm))
treated$tm <- pmin(6, treated$tm)

treated$cenc <- as.numeric(!(treated$centm < treated$tmc))
treated$cenc
treated$tmc[which(treated$cenc == 0)] <- treated$centm[which(treated$cenc == 0)]
treated$cenc <- as.numeric(!(6 < treated$tmc))
treated$tmc <- pmin(6, treated$tmc)





#############
#### Load development PSC (proximity weights + weighted pscfit)

source("/Users/richardjackson/Documents/GitHub/pscDevelop/ProximityWeights/methods/load_psc_dev.R")
#############

#' Summarise simulation output stored in res.array (estimate + SE per replicate).
#' Interval limits use Wald bounds est +/- z * SE (not HPD).
summarise_res_array <- function(res, true_val, method_names = NULL, z = qnorm(0.975)) {
  nr <- dim(res)[1]
  nsim <- dim(res)[3]
  est <- matrix(res[, 1, ], nrow = nr, ncol = nsim)
  se <- matrix(res[, 2, ], nrow = nr, ncol = nsim)

  if (is.null(method_names)) {
    method_names <- seq_len(nr)
  } else if (length(method_names) != nr) {
    stop("method_names length must match nrow(res).")
  }

  lo <- est - z * se
  hi <- est + z * se

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





sa_method_names <- c(
  "Pooled (unadj)", "Pooled (full adj)", "Pooled (partial adj)",
  "PSC full", "PSC partial",
  "SC full", "SC partial",
  "Bayes hist (info)", "Bayes hist (case-weighted)",
  "PSC full (pw)", "PSC partial (pw)"
)

hybrid_method_names <- c(
  "RCT (unadj)", "RCT (adj)",
  "Pooled (unadj)", "Pooled (adj)",
  "PSC full (combined)", "PSC partial (combined)",
  "SC full", "SC partial",
  "PSC full (pw, combined)", "PSC partial (pw, combined)",
  "Bayes vague", "Bayes informative", "Bayes commensurate",
  "Bayes case-weighted", "Bayes partial case-weighted"
)


####################################################################################
####################################################################################
########################## Single Arm Trials #######################################
####################################################################################
####################################################################################

nsim <- 500

res.array <- array(NA, dim = c(11, 2, nsim))
res.array

####
Ncont <- 300
Ntrt <- 50
ns <-1
sim_fail_sa <- rep(NA_character_, nsim)

for (ns in 1:nsim) {

  tryCatch({

  ### Sample
  contDat <- cont[sample(1:nrow(cont), Ncont, replace = F), ]
  trtDat <- treated[sample(1:nrow(treated), Ntrt, replace = F), ]

  ### Defining Time
  contDat$time <- contDat$tmc
  contDat$status <- contDat$cenc
  trtDat$time <- trtDat$tm
  trtDat$status <- trtDat$cen

  ######### Defining Treatment Effects
  contDat$trt <- 0
  trtDat$trt <- 1

  ### Creating combined dataset
  combDat <- rbind(contDat, trtDat)

  #####
  contDat$s.ob <- Surv(contDat$time, contDat$status)
  trtDat$s.ob <- Surv(trtDat$time, trtDat$cen)
  combDat$s.ob <- Surv(combDat$time, combDat$cen)

  ###############
  ### Pooled Analysis
  ###############

  cm <- coxph(s.ob ~ trt, data = combDat)
  cm_adj <- coxph(s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5 + trt, data = combDat)
  cm_padj <- coxph(s.ob ~ X_1 + X_3 + trt, data = combDat)

  us <- survreg(s.ob ~ trt, data = combDat, dist = "exponential")

  ###############
  ### PSC
  ###############

  #### Model
  contDat$X_3 <- factor(contDat$X_3)
  contDat$X_4 <- factor(contDat$X_4)
  contDat$X_5 <- factor(contDat$X_5)

  cfm <- flexsurvspline(
    s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5,
    data = contDat,
    k = 3
  )

  pcfm <- flexsurvspline(
    s.ob ~ X_1 + X_3,
    data = contDat,
    k = 3
  )

  ### development version for proximity weights
  cfmw <- as_pscCFM_dev(cfm)
  pcfmw <- as_pscCFM_dev(pcfm)

  ### Align treated cohort covariates with control CFM factor levels
  trtDat$X_3 <- factor(trtDat$X_3, levels = levels(contDat$X_3))
  trtDat$X_4 <- factor(trtDat$X_4, levels = levels(contDat$X_4))
  trtDat$X_5 <- factor(trtDat$X_5, levels = levels(contDat$X_5))

  ### pscData expects time and cen (not status alone)
  trtDat$cen <- trtDat$status

  ### Estimating weights on full treated cohort (single-arm DC)
  pw <- proximityWeights(cfmw, trtDat, method = "joint")
  ppw <- proximityWeights(pcfmw, trtDat, method = "joint")
  # If method = "all", pass one element to pscfit, e.g. proximity_weights = pw$joint

  #### PSC (unweighted)
  psc_full <- pscfit(cfmw, trtDat, nchain = 1)
  psc_part <- pscfit(pcfmw, trtDat, nchain = 1)

  summary(pw$weights)
  pw$ess

  psc_full$postEst[, c("median", "sd")]
  wpsc_full$postEst[, c("median", "sd")]
  psc_full$postFit[, c("ess_bulk", "rhat")]
  wpsc_full$postFit[, c("ess_bulk", "rhat")]


  #### PSC (proximity-weighted) — pass full pw object, not pw$weights
  wpsc_full <- pscfit(cfmw, trtDat, proximity_weights = pw, nchain = 1)
  wpsc_part <- pscfit(pcfmw, trtDat, proximity_weights = ppw, nchain = 1)

  psc_full
  wpsc_full

  #################
  ### Synthetic Controls
  #################

  X <- as.matrix(combDat[, c("X_1", "X_2", "X_3", "X_4")])
  rm(eb)
  eb <- ebalance(Treatment = combDat$trt, X = X)

  combDat$w <- 1
  combDat$w[combDat$trt == 0] <- eb$w

  pX <- as.matrix(combDat[, c("X_1", "X_3")])
  eb <- ebalance(Treatment = combDat$trt, X = pX)
  combDat$pw <- 1
  combDat$pw[combDat$trt == 0] <- eb$w

  sc_cm <- coxph(s.ob ~ trt, data = combDat, weights = combDat$w)
  sc_pcm <- coxph(s.ob ~ trt, data = combDat, weights = combDat$pw)


  ####################################
  ### Bayesian Historical Controls

  sr <- survreg(s.ob ~ 1, data = contDat, dist = "exponential")

  co <- as.numeric(sr$coefficients)
  se <- as.numeric(sqrt(sr$var))

  stan_data <- list(
    N = nrow(trtDat),
    time = trtDat$time,
    status = trtDat$cen,
    x = trtDat$trt,
    lam0_mn = co,
    lam0_t = se
  )

  singArm_baye <- stan(
    model_code = exp_code,
    data = stan_data,
    iter = 2000,
    chains = 2,
    seed = 21319,
    refresh = 0
  )

  ######## Case Weighted
  cp_w <- combDat$w[combDat$trt == 0]

  stan_cw_data <- list(
    N = nrow(trtDat),
    time = trtDat$time,
    status = trtDat$cen,
    x = trtDat$trt,
    N0 = nrow(contDat),
    timec = contDat$time,
    statusc = contDat$cen,
    a0 = cp_w
  )

  singArm_caseW_baye <- stan(
    model_code = exp_caseW_code,
    data = stan_cw_data,
    iter = 2000,
    chains = 2,
    seed = 21319,
    refresh = 0
  )

  ### Collecting results (col 1 = log(HR) estimate, col 2 = SE)
  r1 <- summary(cm)$coef[1, c(1, 3)]
  r2 <- summary(cm_adj)$coef[6, c(1, 3)]
  r3 <- summary(cm_padj)$coef[3, c(1, 3)]

  r4 <- c(as.numeric(psc_full$postEst$median[1]), as.numeric(psc_full$postEst$sd[1]))
  r5 <- c(as.numeric(psc_part$postEst$median[1]), as.numeric(psc_part$postEst$sd[1]))

  r6 <- summary(sc_cm)$coef[1, c(1, 3)]
  r7 <- summary(sc_pcm)$coef[1, c(1, 3)]

  r8 <- summary(singArm_baye)$summary[2, c(1, 3)]
  r9 <- summary(singArm_caseW_baye)$summary[2, c(1, 3)]

  r10 <- c(as.numeric(wpsc_full$postEst$median[1]), as.numeric(wpsc_full$postEst$sd[1]))
  r11 <- c(as.numeric(wpsc_part$postEst$median[1]), as.numeric(wpsc_part$postEst$sd[1]))

  tmp.res <- rbind(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11)
  res.array[,, ns] <- tmp.res

  plot(rep(c(1:11), nsim), c(res.array[, 1, ]), ylim = c(-2, 2), pch = 20, col = c(1:11), main = ns,
       xlab="log(HR)",ylab="scenario")

  }, error = function(e) {
    sim_fail_sa[ns] <<- conditionMessage(e)
    warning("Single-arm iteration ", ns, " failed: ", conditionMessage(e), call. = FALSE)
  })

}

cat(
  "Single-arm: ", sum(!is.na(sim_fail_sa)), " of ", nsim, " iterations failed\n",
  sep = ""
)


sim_summary_sa <- summarise_res_array(res.array, true_val = beta, method_names = sa_method_names)
sim_summary_sa

setwd("~/Documents/GitHub/synthetic-Controlled-Trials/Examples/Data/SimRes_td")

fn <- paste("res.array_sa_c", Ncont, "_e", Ntrt, ".R", sep = "")
save(res.array, sim_fail_sa, sim_summary_sa, file = fn)








################################################################################################
################################################################################################
######################################## Hybrid Trials #########################################
################################################################################################
################################################################################################

##### Hybrid trial (1:1 randomisation)

nsim <- 50

rand.res.array <- array(NA, dim = c(15, 2, nsim))
rand.res.array

sim_fail_hybrid <- rep(NA_character_, nsim)

for (ns in 1:nsim) {

  tryCatch({

  ### Sample
  contDat <- cont[sample(1:nrow(cont), Ncont, replace = F), ]
  trtDat <- treated[sample(1:nrow(treated), Ntrt, replace = F), ]

  ######## Defining treatment
  contDat$trt <- 0
  trtDat$trt <- rbinom(Ntrt, 1, 0.5)

  ### Defining time and censoring (align with single-arm setup)
  contDat$time <- contDat$tmc
  contDat$status <- contDat$cenc

  trtDat$time <- trtDat$tm
  trtDat$status <- trtDat$cen
  trtDat$time[trtDat$trt == 0] <- trtDat$tmc[trtDat$trt == 0]
  trtDat$status[trtDat$trt == 0] <- trtDat$cenc[trtDat$trt == 0]

  ### Creating combined dataset
  contDat$source <- "hist"
  trtDat$source <- "cont"
  combDat <- rbind(contDat, trtDat)
  combDat$cen <- ifelse(combDat$source == "hist", combDat$cenc, combDat$cen)

  contDat$s.ob <- Surv(contDat$time, contDat$status)
  trtDat$s.ob <- Surv(trtDat$time, trtDat$status)
  combDat$s.ob <- Surv(combDat$time, combDat$cen)

  #### Standard randomised comparison
  cm.rand <- coxph(s.ob ~ trt, data = trtDat)
  cm.rand.ad <- coxph(s.ob ~ X_1 + X_3 + trt, data = trtDat)

  #### Pooled analysis
  cm.pooled <- coxph(s.ob ~ trt, data = combDat)
  cm.pooled.ad <- coxph(s.ob ~ X_1 + X_3 + trt, data = combDat)

  ###############
  ### PSC
  ###############

  contDat$X_3 <- factor(contDat$X_3)
  contDat$X_4 <- factor(contDat$X_4)
  contDat$X_5 <- factor(contDat$X_5)

  cfm <- flexsurvspline(
    s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5,
    data = contDat,
    k = 3
  )

  pcfm <- flexsurvspline(
    s.ob ~ X_1 + X_3,
    data = contDat,
    k = 3
  )

  ### development version for proximity weights
  cfmw <- as_pscCFM_dev(cfm)
  pcfmw <- as_pscCFM_dev(pcfm)

  ### Align trial cohort covariates with control CFM factor levels
  trtDat$X_3 <- factor(trtDat$X_3, levels = levels(contDat$X_3))
  trtDat$X_4 <- factor(trtDat$X_4, levels = levels(contDat$X_4))
  trtDat$X_5 <- factor(trtDat$X_5, levels = levels(contDat$X_5))

  ### pscData expects time and cen (not status alone)
  trtDat$cen <- trtDat$status

  trtDat$trt <- factor(trtDat$trt + 1)

  ### Estimating weights on full trial cohort (hybrid DC)
  pw <- proximityWeights(cfmw, trtDat, method = "joint")
  ppw <- proximityWeights(pcfmw, trtDat, method = "joint")
  # If method = "all", pass one element to pscfit, e.g. proximity_weights = pw$joint

  #### PSC (unweighted)
  psc_full <- pscfit(cfmw, trtDat, nchain = 1, trt = trtDat$trt)
  psc_part <- pscfit(pcfmw, trtDat, nchain = 1, trt = trtDat$trt)

  #### PSC (proximity-weighted) — pass full pw object, not pw$weights
  wpsc_full <- pscfit(cfmw, trtDat, proximity_weights = pw, nchain = 1, trt = trtDat$trt)
  wpsc_part <- pscfit(pcfmw, trtDat, proximity_weights = ppw, nchain = 1, trt = trtDat$trt)

  psc_fullc <- pscComb(psc_full)
  psc_partc <- pscComb(psc_part)
  wpsc_fullc <- pscComb(wpsc_full)
  wpsc_partc <- pscComb(wpsc_part)

  #################
  ### Synthetic controls
  #################
  contCont <- combDat[which(combDat$trt == 0 & combDat$source == "cont"), ]
  synthDat <- combDat[-which(combDat$trt == 0 & combDat$source == "cont"), ]

  X <- as.matrix(synthDat[, c("X_1", "X_2", "X_3", "X_4", "X_5")])
  eb <- ebalance(Treatment = synthDat$trt, X = X)

  pX <- as.matrix(synthDat[, c("X_1", "X_3")])
  peb <- ebalance(Treatment = synthDat$trt, X = pX)

  synthDat$w <- 1
  synthDat$w[synthDat$trt == 0] <- eb$w

  synthDat$pw <- 1
  synthDat$pw[synthDat$trt == 0] <- peb$w
  contCont$w <- 1
  contCont$pw <- 1

  scData <- rbind(synthDat, contCont)
  scData$s.ob <- Surv(scData$time, scData$cen)
  sc_cm <- coxph(s.ob ~ trt, data = scData, weights = scData$w)
  sc_pcm <- coxph(s.ob ~ trt, data = scData, weights = scData$pw)

  ####################################
  ### Bayesian historical controls
  ####################################

  sr <- survreg(s.ob ~ 1, data = contDat, dist = "exponential")
  co <- as.numeric(sr$coefficients)
  se <- as.numeric(sqrt(sr$var))

  x <- as.numeric(trtDat$trt) - 1

  vague_data <- list(
    N = nrow(trtDat),
    time = trtDat$time,
    status = trtDat$status,
    x = x,
    lam0_mn = 0,
    lam0_t = 10
  )

  vague_baye <- stan(
    model_code = exp_code,
    data = vague_data,
    iter = 2000,
    chains = 2,
    seed = 21319,
    refresh = 0
  )

  stan_data <- list(
    N = nrow(trtDat),
    time = trtDat$time,
    status = trtDat$status,
    x = x,
    lam0_mn = co,
    lam0_t = se
  )

  hybrid_baye <- stan(
    model_code = exp_code,
    data = stan_data,
    iter = 2000,
    chains = 2,
    seed = 21319,
    refresh = 0
  )

  stan_comm_data <- list(
    N = nrow(trtDat),
    time = trtDat$time,
    status = trtDat$status,
    x = x,
    lam0_mn = co
  )

  hybrid_cp_baye <- stan(
    model_code = exp_comm_code,
    data = stan_comm_data,
    iter = 2000,
    chains = 2,
    seed = 21319,
    refresh = 0
  )

  CW_Data <- synthDat[which(synthDat$source == "hist"), ]
  cw_w <- CW_Data$w

  stan_cw_data <- list(
    N = nrow(trtDat),
    time = trtDat$time,
    status = trtDat$status,
    x = x,
    N0 = nrow(CW_Data),
    timec = CW_Data$time,
    statusc = CW_Data$cen,
    a0 = cw_w
  )

  hybrid_caseW_baye <- stan(
    model_code = exp_caseW_code,
    data = stan_cw_data,
    iter = 2000,
    chains = 2,
    seed = 21319,
    refresh = 0
  )

  cw_pw <- CW_Data$pw
  stan_cw_data <- list(
    N = nrow(trtDat),
    time = trtDat$time,
    status = trtDat$status,
    x = x,
    N0 = nrow(contDat),
    timec = contDat$time,
    statusc = contDat$status,
    a0 = cw_pw
  )

  hybrid_pcaseW_baye <- stan(
    model_code = exp_caseW_code,
    data = stan_cw_data,
    iter = 2000,
    chains = 2,
    seed = 21319,
    refresh = 0
  )

  ### Collecting results (col 1 = log(HR) estimate, col 2 = SE)
  r1 <- summary(cm.rand)$coef[1, c(1, 3)]
  r2 <- summary(cm.rand.ad)$coef[3, c(1, 3)]

  r3 <- summary(cm.pooled)$coef[1, c(1, 3)]
  r4 <- summary(cm.pooled.ad)$coef[3, c(1, 3)]

  r5 <- c(median(as.numeric(unlist(psc_fullc))), sd(as.numeric(unlist(psc_fullc))))
  r6 <- c(median(as.numeric(unlist(psc_partc))), sd(as.numeric(unlist(psc_partc))))

  r7 <- summary(sc_cm)$coef[1, c(1, 3)]
  r8 <- summary(sc_pcm)$coef[1, c(1, 3)]

  r9 <- c(median(as.numeric(unlist(wpsc_fullc))), sd(as.numeric(unlist(wpsc_fullc))))
  r10 <- c(median(as.numeric(unlist(wpsc_partc))), sd(as.numeric(unlist(wpsc_partc))))

  r11 <- summary(vague_baye)$summary[2, c(1, 3)]
  r12 <- summary(hybrid_baye)$summary[2, c(1, 3)]
  r13 <- summary(hybrid_cp_baye)$summary[2, c(1, 3)]
  r14 <- summary(hybrid_caseW_baye)$summary[2, c(1, 3)]
  r15 <- summary(hybrid_pcaseW_baye)$summary[2, c(1, 3)]

  tmp.res <- rbind(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15)
  rand.res.array[,, ns] <- tmp.res
  rand.res.array

  plot(rep(c(1:15), nsim), c(rand.res.array[, 1, ]), ylim = c(-2, 2), pch = 20, col = c(1:15), main = ns,
       xlab = "log(HR)", ylab = "scenario")

  }, error = function(e) {
    sim_fail_hybrid[ns] <<- conditionMessage(e)
    warning("Hybrid iteration ", ns, " failed: ", conditionMessage(e), call. = FALSE)
  })

}

cat(
  "Hybrid: ", sum(!is.na(sim_fail_hybrid)), " of ", nsim, " iterations failed\n",
  sep = ""
)

sim_summary_hybrid <- summarise_res_array(
  rand.res.array,
  true_val = beta,
  method_names = hybrid_method_names
)
print_res_summary(sim_summary_hybrid, "Hybrid trial simulation summary", beta)

setwd("~/Documents/GitHub/synthetic-Controlled-Trials/Examples/Data/SimRes_td")

fn <- paste("res.array_ra_c", Ncont, "_e", Ntrt, ".R", sep = "")
save(rand.res.array, sim_fail_hybrid, sim_summary_hybrid, file = fn)


################################################################################################
################################################################################################
########################################## Stan Code ###########################################
################################################################################################
################################################################################################

exp_code <- "
data {
  int<lower=0> N;
  vector<lower=0>[N] time;
  real lam0_mn;
  real lam0_t;
  int<lower=0,upper=1> status[N];
  vector[N] x;  // covariate: treatment
}

parameters {
  real<lower=0> lambda0;
  real beta;
}

model {
  // Priors
  lambda0 ~ normal(lam0_mn,lam0_t);
  beta ~ normal(0,10);

  // Likelihood with right censoring
  for (i in 1:N) {
    real lambda_i;
    lambda_i =  exp(-lambda0+beta * x[i]);

    if (status[i] == 1)
      target += exponential_lpdf(time[i] | lambda_i);
    else
      target += exponential_lccdf(time[i] | lambda_i);
  }
}"


exp_comm_code <- "
data {

// --- Current RCT Data ---
  int<lower=0> N;
  vector<lower=0>[N] time;
  int<lower=0,upper=1> status[N];
  vector[N] x;
  real lam0_mn;

}

parameters {
  real lambda0;           // baseline log-hazard
  real beta;            // treatment effect
  real<lower=0> tau;    // commensurability (precision)
}

model {
  // Priors
  lambda0 ~ normal(lam0_mn, 1/sqrt(tau));
  beta ~ normal(0, 10);
  tau ~ gamma(1, 1);

  // Likelihood
  for (i in 1:N) {
    real lambda_i;
    lambda_i = exp(-lambda0 + beta * x[i]);

    if (status[i] == 1)
      target += exponential_lpdf(time[i] | lambda_i);
    else
      target += exponential_lccdf(time[i] | lambda_i);
  }
}
"


exp_caseW_code <- "
data {
  int<lower=0> N;
  vector<lower=0>[N] time;
  int<lower=0,upper=1> status[N];
  vector[N] x;

  // --- RWD Data ---
  int<lower=0> N0;
  vector<lower=0>[N0] timec;
  int<lower=0,upper=1> statusc[N0];
  vector[N0] a0;


}

parameters {
  real lambda0;           // baseline log-hazard
  real beta;            // treatment effect
}

model {
  // Priors
  lambda0 ~ normal(0, 10);
  beta ~ normal(0, 10);


  // RCT Likelihood
  for (i in 1:N) {
    real lambda_i;
    lambda_i = exp(-lambda0 + beta * x[i]);

    if (status[i] == 1)
      target += exponential_lpdf(time[i] | lambda_i);
    else
      target += exponential_lccdf(time[i] | lambda_i);
  }

  // RWD Likelihood
   for (j in 1:N0) {
    real lambdac_j;
    lambdac_j = exp(-lambda0);

    if (statusc[j] == 1)
      target += a0[j] * exponential_lpdf(timec[j] | lambdac_j);
    else
      target += a0[j] * exponential_lccdf(timec[j] | lambdac_j);
  }
}
"


