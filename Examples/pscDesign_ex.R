### Comparing, SC, PSC and pooled analysis

###
library(survival)
library(psc)
library(flexsurv)


######## CONTROL POPULATION ############
### Creating population of control and experimental patients
N <- 5000

gam_0 <- 0.05
gam_1 <- log(0.6)
gam_2 <- log(0.7)
gam_3 <- log(0.6)
gam_4 <- log(0.7)
gam_5 <- log(0.8)
beta <- log(0.7)

X_1 <- rnorm(N)
X_2 <- rnorm(N)
X_3 <- rbinom(N,1,0.3)
X_4 <- rbinom(N,1,0.4)
X_5 <- rbinom(N,1,0.5)

trt <- rbinom(N,1,0.5)

lp <- exp(gam_0+gam_1*X_1 + gam_2*X_2 + gam_3*X_3 +
            gam_4*X_4 + gam_5*X_5 + beta)
lpc <- exp(gam_0+gam_1*X_1 + gam_2*X_2 + gam_3*X_3 +
             gam_4*X_4 + gam_5*X_5)

tm <- rexp(N,lp)
tmc <- rexp(N,lpc)
centm <- rexp(N,rep(0.02,N));centm

### Creating dataset
cont <- data.frame(X_1,X_2,X_3,X_4,X_5,tm,tmc,centm)


######## Treated POPULATION ############
### Creating population of control and experimental patients
N <- 5000

gam_0 <- 0.05
gam_1 <- log(0.6)
gam_2 <- log(0.7)
gam_3 <- log(0.6)
gam_4 <- log(0.7)
gam_5 <- log(0.8)
beta <- log(0.7)

X_1 <- rnorm(N)*1.25
X_2 <- rnorm(N)*0.75
X_3 <- rbinom(N,1,0.5)
X_4 <- rbinom(N,1,0.6)
X_5 <- rbinom(N,1,0.7)

lp <- exp(gam_0+gam_1*X_1 + gam_2*X_2 + gam_3*X_3 +
            gam_4*X_4 + gam_5*X_5 + beta)
lpc <- exp(gam_0+gam_1*X_1 + gam_2*X_2 + gam_3*X_3 +
            gam_4*X_4 + gam_5*X_5)

tm <- rexp(N,lp)
tmc <- rexp(N,lpc)
centm <- rexp(N,rep(0.02,N));centm

### Creating dataset
treated <- data.frame(X_1,X_2,X_3,X_4,X_5,tm,tmc,centm)



#### Consoring patterns
cont$cen <- as.numeric(!(cont$centm<cont$tm))
cont$tm[which(cont$cen==0)] <- cont$centm[which(cont$cen==0)]
cont$cenc <- as.numeric(!(cont$centm<cont$tmc))
cont$tmc[which(cont$cenc==0)] <- cont$centm[which(cont$cenc==0)]

treated$cen <- as.numeric(!(treated$centm<treated$tm));treated$cen
treated$tm[which(treated$cen==0)] <- treated$centm[which(treated$cen==0)]
treated$cen <- as.numeric(!(6<treated$tm))
treated$tm <- pmin(6,treated$tm)

treated$cenc <- as.numeric(!(treated$centm<treated$tmc));treated$cenc
treated$tmc[which(treated$cenc==0)] <- treated$centm[which(treated$cenc==0)]
treated$cenc <- as.numeric(!(6<treated$tmc))
treated$tmc <- pmin(6,treated$tmc)



############# Single Arm Trials

####
Ncont <- 500
Ntrt <- 100

### Sample
contDat <- cont[sample(1:nrow(cont),Ncont,replace=F),]
trtDat<- treated[sample(1:nrow(treated),Ntrt,replace=F),]

#####
contDat$s.ob <- Surv(contDat$tm,contDat$cen)
trtDat$s.ob <- Surv(trtDat$tm,trtDat$cen)

######### Pooled Analysis
contDat$trt <- 0
trtDat$trt <- 1

combData <- rbind(contDat,trtDat)
combData$s.ob <- Surv(combData$tm,combData$cen)

cm <- coxph(s.ob~trt,data=combData)
cm_adj <- coxph(s.ob~X_1+X_2+X_3+X_4+X_5+trt,data=combData)
cm_padj <- coxph(s.ob~X_1+X_3+trt,data=combData)

us <- survreg(s.ob~trt,data=combData,dist="exponential")
summary(us)
summary(cm)



### PSC

#### Model
contDat$time <- contDat$tm
contDat$s.ob <- Surv(contDat$time,contDat$cen)

cfm <- flexsurvspline(s.ob~X_1+X_2+X_3+X_4+X_5,data=contDat,k=3)
pcfm <- flexsurvspline(s.ob~X_1+X_3,data=contDat,k=3)


trtDat$time <- trtDat$tm

psc_full <- pscfit(cfm,trtDat)
psc_part <- pscfit(pcfm,trtDat)


#### Bayesian Analysis
library(rstan)

sr <- survreg(s.ob~1,data=contDat,dist="exponential");sr

co <- as.numeric(sr$coefficients)
se <- as.numeric(sqrt(sr$var))


stan_data <- list(N=nrow(trtDat),time=trtDat$time,status=trtDat$cen,
                  x=trtDat$trt,lam0_mn=co,lam0_t=se)

fit <- stan(model_code=exp_code,data=stan_data,  iter = 2000,chains = 2,
            seed = 21319,refresh=0)




############################
##### Hybrid Trial (1:1) randomisation

####
Ncont <- 250
Ntrt <- 100



### Sample
contDat <- cont[sample(1:nrow(cont),Ncont,replace=F),]
trtDat<- treated[sample(1:nrow(treated),Ntrt,replace=F),]

##
trtDat$trt <- rbinom(Ntrt,1,0.66)
trtDat$time <- trtDat$tm
trtDat$time[trtDat$trt==0] <- trtDat$tmc[trtDat$trt==0]

cm <-coxph(Surv(time,cen)~trt,data=trtDat)


#### Model
contDat$time <- contDat$tm
contDat$s.ob <- Surv(contDat$time,contDat$cen)

cfm <- flexsurvspline(s.ob~X_1+X_2+X_3+X_4+X_5,data=contDat,k=3)
pcfm <- flexsurvspline(s.ob~X_1+X_3,data=contDat,k=3)

id <- which(trtDat$trt==1)

psc_full <- pscfit(cfm,trtDat[id,])
psc_part <- pscfit(pcfm,trtDat[id,])
plot(psc_full)
plot(psc_part)

### Combining
cm.co <- summary(cm)$coefficients[c(1,3)]
psc.co <- coef(psc_full)

post.psc.full.sd <- 1/(1/as.numeric(psc.co[3]) + 1/as.numeric(cm.co[2]))
post.psc.full.mu <- ((as.numeric(psc.co[2])/as.numeric(psc.co[3])) + (as.numeric(cm.co[1])/as.numeric(cm.co[2])))*post.sd

psc.co <- coef(psc_part)
post.psc.part.sd <- 1/(1/as.numeric(psc.co[3]) + 1/as.numeric(cm.co[2]))
post.psc.part.mu <- ((as.numeric(psc.co[2])/as.numeric(psc.co[3])) + (as.numeric(cm.co[1])/as.numeric(cm.co[2])))*post.sd




#### Commensurate Prior


#### Bayesian Analysis

sr <- survreg(s.ob~1,data=contDat,dist="exponential");sr

co <- as.numeric(sr$coefficients);co
se <- as.numeric(sqrt(sr$var));se

stan_data <- list(N=nrow(trtDat),time=trtDat$time,status=trtDat$cen,
                  x=trtDat$trt,lam0_mn=co,lam0_t=se)

fit_baye <- stan(model_code=exp_code,data=stan_data,  iter = 2000,chains = 2,
            seed = 21319,refresh=0)


fit_baye

stan_comm_data <- list(N=nrow(trtDat),time=trtDat$time,status=trtDat$cen,
                       x=trtDat$trt,lam0_mn=co)

fit_cp <- stan(model_code=exp_comm_code,data=stan_data,  iter = 2000,chains = 2,
               seed = 21319,refresh=0)


cm
fit_baye
fit_cp
post.psc.full.mu
post.psc.full.sd
post.psc.part.mu
post.psc.part.sd
psc_full
cm

psc_full



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
  int<lower=0> N;
  vector<lower=0>[N] time;
  int<lower=0,upper=1> status[N];
  real lam0_mn;
  vector[N] x;
}

parameters {
  real lambda0;           // baseline log-hazard
  real beta;            // treatment effect
  real<lower=0> tau;    // commensurability (precision)
}

model {
  // Priors
  lambda0 ~ normal(lam0_mn,1 / sqrt(tau));
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










