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







####################################################################################
####################################################################################
########################## Single Arm Trials #######################################
####################################################################################
####################################################################################

nsim <- 1000

res.array <- array(NA, dim = c(9, 2, nsim))
res.array


ns <- 1

for (ns in 1:nsim) {

    ####
  Ncont <- 500
  Ntrt <- 100

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

  combDat[1:3,]

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


  psc_full <- pscfit(cfm, trtDat, nchain = 1)
  psc_part <- pscfit(pcfm, trtDat, nchain = 1)

  #################
  ### Synthetic Controls
  #################

  X <- as.matrix(combDat[, c("X_1", "X_2", "X_3", "X_4", "X_5")])
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

  ### Collecting Results
  r1 <- summary(cm)$coef[1, c(1, 3)]
  r2 <- summary(cm_adj)$coef[6, c(1, 3)]
  r3 <- summary(cm_padj)$coef[3, c(1, 3)]

  r4 <- as.numeric(coef(psc_full)[c(2, 3)])
  r5 <- as.numeric(coef(psc_part)[c(2, 3)])

  r6 <- summary(sc_cm)$coef[c(1, 3)]
  r7 <- summary(sc_pcm)$coef[c(1, 3)]

  r8 <- summary(singArm_baye)$summary[2, c(1, 3)]
  r9 <- summary(singArm_caseW_baye)$summary[2, c(1, 3)]

  tmp.res <- rbind(r1, r2, r3, r4, r5, r6, r7, r8, r9)
  res.array[,, ns] <- tmp.res
  res.array
  plot(rep(c(1:9),nsim),c(res.array[,1,]),ylim=c(-2,2),pch=20,col=c(1:9),main=ns,
       xlab="log(HR)",ylab="scenario")

}





dim(res.array)


plot(rep(c(1:9),50),c(res.array[,1,]),ylim=c(-2,2),pch=20,col=c(1:9),main=ns,
     xlab="log(HR)",ylab="scenario")
abline(h=log(0.7))

log(0.7)
rowMeans(res.array[,1,],na.rm=T)



res.array[,1]





################################################################################################
################################################################################################
######################################## Hybrid Trials #########################################
################################################################################################
################################################################################################

##### Hybrid Trial (1:1) randomisation

####
Ncont <- 500
Ntrt <- 100


### Sample
contDat <- cont[sample(1:nrow(cont), Ncont, replace = F), ]
trtDat <- treated[sample(1:nrow(treated), Ntrt, replace = F), ]

######## Defining Treatment
contDat$trt <- 0
trtDat$trt <- rbinom(Ntrt, 1, 0.5)


##
contDat$time <- contDat$tmc
trtDat$time <- trtDat$tm
trtDat$time[trtDat$trt == 0] <- trtDat$tmc[trtDat$trt == 0]
contDat$trt <- 0
contDat$time <- contDat$tmc

### Creating combined dataset
contDat$source <- "hist"
trtDat$source <- "cont"
combDat <- rbind(contDat, trtDat)

#####
contDat$s.ob <- Surv(contDat$time, contDat$cen)
trtDat$s.ob <- Surv(trtDat$time, trtDat$cen)
combDat$s.ob <- Surv(combDat$time, combDat$cen)


#### Standard Randomised Comparison
cm.rand <- coxph(s.ob ~ trt, data = trtDat)
cm.rand.ad <- coxph(s.ob ~ X_1 + X_3 + trt, data = trtDat)

trtDat[1:3, ]

#### Pooled Analysis
cm.pooled <- coxph(s.ob ~ trt, data = combDat)


######################################
########### Personalised Synthetic Controls
######################################

### contDat$time <- contDat$tm
contDat$s.ob <- Surv(contDat$time, contDat$cen)

cfm <- flexsurvspline(s.ob ~ X_1 + X_2 + X_3 + X_4 + X_5, data = contDat, k = 3)
pcfm <- flexsurvspline(s.ob ~ X_1 + X_3, data = contDat, k = 3)

tid <- which(trtDat$trt == 2)
trtDat$trt <- factor(trtDat$trt + 1)

psc_full <- pscfit(cfm, trtDat, nchain = 1, trt = trtDat$trt, nsim = 10000)
psc_part <- pscfit(pcfm, trtDat, nchain = 1, trt = trtDat$trt, nsim = 10000)


cm.rand.ad

psc_full

pscComb(psc_full)
pscComb(psc_part)


######################################
######################################

#################
### Synthetic Controls
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


#plot(psc_full)
#plot(psc_part)

#### Commensurate Prior
#### Bayesian Analysis

sr <- survreg(s.ob ~ 1, data = contDat, dist = "exponential")
co <- as.numeric(sr$coefficients)
se <- as.numeric(sqrt(sr$var))


vague_data <- list(
  N = nrow(trtDat),
  time = trtDat$time,
  status = trtDat$cen,
  x = trtDat$trt,
  lam0_mn = 0,
  lam0_t = 10
)

stan_data <- list(
  N = nrow(trtDat),
  time = trtDat$time,
  status = trtDat$cen,
  x = trtDat$trt,
  lam0_mn = co,
  lam0_t = se
)

stan_comm_data <- list(
  N = nrow(trtDat),
  time = trtDat$time,
  status = trtDat$cen,
  x = trtDat$trt,
  lam0_mn = co
)


vague_baye <- stan(
  model_code = exp_code,
  data = vague_data,
  iter = 2000,
  chains = 2,
  seed = 21319,
  refresh = 0
)

hybrid_baye <- stan(
  model_code = exp_code,
  data = stan_data,
  iter = 2000,
  chains = 2,
  seed = 21319,
  refresh = 0
)

hybrid_cp_baye <- stan(
  model_code = exp_comm_code,
  data = stan_comm_data,
  iter = 2000,
  chains = 2,
  seed = 21319,
  refresh = 0
)


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


hybrid_caseW_baye <- stan(
  model_code = exp_caseW_code,
  data = stan_cw_data,
  iter = 2000,
  chains = 2,
  seed = 21319,
  refresh = 0
)


cm.rand
sc_cm
sc_pcm

pscComb(psc_full)
pscComb(psc_part)


hybrid_cp_baye

exp_caseW_code


cm.pooled
cm.rand
sc_cm
vague_baye
co
hybrid_baye
co
hybrid_cp_baye
hybrid_caseW_baye


cp_w <- synthDat$w[synthDat$trt == 0]


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


