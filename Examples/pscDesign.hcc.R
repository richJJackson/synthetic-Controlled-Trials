

### HCC Example

#install.packages("pscDesign")

setwd("/Users/richardjackson/Documents/GitHub/pscDesign")
devtools::load_all()

setwd("/Users/richardjackson/Documents/GitHub/psc")
devtools::load_all()

#library(psc)
#library(pscDesign)


### Getting Atezo Bev model

setwd("/Users/richardjackson/Documents/GitHub/synthetic-Controlled-Trials/Examples/Data")

load("abCFM_prog.Rda")


### simulation paramteres (Once all scenarios complete can you increase these to 500, 2000, 400)
nsim = 200
nsim.psc = 1000
burn.psc = 200

### Setting 'delta'
beta <- log(0.7)

### Setting recruitment
nSite <- 10
openRate <- 1
maxTime <- 18



#######################################
### 6 different recruitment stratergies based on 40 - 150 patients

rpm1 <- 0.31
rec1 <- recForcast(nSite,rpm1,openRate,Max.Time = maxTime)

rpm2 <- 0.385
rec2 <- recForcast(nSite,rpm2,openRate,Max.Time = maxTime)

rpm3 <- 0.5
rec3 <- recForcast(nSite,rpm3,openRate,Max.Time = maxTime)

rpm4 <- 0.58
rec4 <- recForcast(nSite,rpm4,openRate,Max.Time = maxTime)

rpm5 <- 0.77
rec5 <- recForcast(nSite,rpm5,openRate,Max.Time = maxTime)

rpm6 <- 1.15
rec6 <- recForcast(nSite,rpm6,openRate,Max.Time = maxTime)

#######################################




#############################

#######################################
### single Arm design
design_singArm1 <- pscDesign(abCFM_prog,beta=beta,n0=0,n1=40,fuTime=12,rec=rec1,nsim=nsim)
design_singArm2 <- pscDesign(abCFM_prog,beta=beta,n0=0,n1=50,fuTime=12,rec=rec2,nsim=nsim)
design_singArm3 <- pscDesign(abCFM_prog,beta=beta,n0=0,n1=65,fuTime=12,rec=rec3,nsim=nsim)
design_singArm4 <- pscDesign(abCFM_prog,beta=beta,n0=0,n1=75,fuTime=12,rec=rec4,nsim=nsim)
design_singArm5 <- pscDesign(abCFM_prog,beta=beta,n0=0,n1=100,fuTime=12,rec=rec5,nsim=nsim)
design_singArm6 <- pscDesign(abCFM_prog,beta=beta,n0=0,n1=150,fuTime=12,rec=rec6,nsim=nsim)

#######################################
### Randomised (1:1)
design_randEq1 <- pscDesign(abCFM_prog,beta=beta,n0=20,n1=20,fuTime=12,rec=rec1,nsim=nsim)
design_randEq2 <- pscDesign(abCFM_prog,beta=beta,n0=25,n1=25,fuTime=12,rec=rec2,nsim=nsim)
design_randEq3 <- pscDesign(abCFM_prog,beta=beta,n0=32,n1=33,fuTime=12,rec=rec3,nsim=nsim)
design_randEq4 <- pscDesign(abCFM_prog,beta=beta,n0=37,n1=38,fuTime=12,rec=rec4,nsim=nsim)
design_randEq5 <- pscDesign(abCFM_prog,beta=beta,n0=50,n1=50,fuTime=12,rec=rec5,nsim=nsim)
design_randEq6 <- pscDesign(abCFM_prog,beta=beta,n0=75,n1=150,fuTime=12,rec=rec6,nsim=nsim)

### Randomised (1:2)
design_randUneq1 <- pscDesign(abCFM_prog,beta=beta,n0=13,n1=27,fuTime=12,rec=rec1,nsim=nsim)
design_randUneq2 <- pscDesign(abCFM_prog,beta=beta,n0=17,n1=33,fuTime=12,rec=rec2,nsim=nsim)
design_randUneq3 <- pscDesign(abCFM_prog,beta=beta,n0=22,n1=43,fuTime=12,rec=rec3,nsim=nsim)
design_randUneq4 <- pscDesign(abCFM_prog,beta=beta,n0=25,n1=50,fuTime=12,rec=rec4,nsim=nsim)
design_randUneq5 <- pscDesign(abCFM_prog,beta=beta,n0=33,n1=67,fuTime=12,rec=rec5,nsim=nsim)
design_randUneq6 <- pscDesign(abCFM_prog,beta=beta,n0=50,n1=100,fuTime=12,rec=rec6,nsim=nsim)

### Randomised (1:3)
design_rand3Uneq1 <- pscDesign(abCFM_prog,beta=beta,n0=10,n1=30,fuTime=12,rec=rec1,nsim=nsim)
design_rand3Uneq2 <- pscDesign(abCFM_prog,beta=beta,n0=12,n1=38,fuTime=12,rec=rec2,nsim=nsim)
design_rand3Uneq3 <- pscDesign(abCFM_prog,beta=beta,n0=16,n1=49,fuTime=12,rec=rec3,nsim=nsim)
design_rand3Uneq4 <- pscDesign(abCFM_prog,beta=beta,n0=19,n1=56,fuTime=12,rec=rec4,nsim=nsim)
design_rand3Uneq5 <- pscDesign(abCFM_prog,beta=beta,n0=25,n1=75,fuTime=12,rec=rec5,nsim=nsim)
design_rand3Uneq6 <- pscDesign(abCFM_prog,beta=beta,n0=30,n1=120,fuTime=12,rec=rec6,nsim=nsim)







#############################. Summarising results

res_sing <- cbind(design_singArm1[[2]][,2],design_singArm2[[2]][,2],design_singArm3[[2]][,2],
             design_singArm4[[2]][,2],design_singArm5[[2]][,2],design_singArm6[[2]][,2])


res_randEq <- cbind(design_randEq1[[2]][,2],design_randEq2[[2]][,2],design_randEq3[[2]][,2],
                  design_randEq4[[2]][,2],design_randEq5[[2]][,2],design_randEq6[[2]][,2])

res_randUneq <- cbind(design_randUneq1[[2]][,2],design_randUneq2[[2]][,2],design_randUneq3[[2]][,2],
                   design_randUneq4[[2]][,2],design_randUneq5[[2]][,2],design_randUneq6[[2]][,2])

res_rand3Uneq <- cbind(design_rand3Uneq1[[2]][,2],design_rand3Uneq2[[2]][,2],design_rand3Uneq3[[2]][,2],
                           design_rand3Uneq4[[2]][,2],design_rand3Uneq5[[2]][,2],design_rand3Uneq6[[2]][,2])

design_rand3Uneq2
res_rand3Uneq

par(mfrow=c(1,3))

plot(ss,res_sing[3,],typ="o",col=1,lwd=2,pch=15,ylim=c(0.3,1),xlim=c(40,150),main="alpha=0.05",
     ylab="Sample Siz",xlab="Power")
lines(ss,res_randEq[3,],typ="o",col=2,lwd=2,pch=15,ylim=c(0.5,1),xlim=c(40,150))
lines(ss,res_randUneq[3,],typ="o",col=3,lwd=2,pch=15,ylim=c(0.5,1),xlim=c(40,150))
lines(ss,res_rand3Uneq[3,],typ="o",col=4,lwd=2,pch=15,ylim=c(0.5,1),xlim=c(40,150))
abline(h=c(0.8,0.9),lty=3,col=6,lwd=2)


plot(ss,res_sing[4,],typ="o",col=1,lwd=2,pch=15,ylim=c(0.3,1),xlim=c(40,150),main="alpha=0.1",
     ylab="Sample Siz",xlab="Power")
lines(ss,res_randEq[4,],typ="o",col=2,lwd=2,pch=15,ylim=c(0.5,1),xlim=c(40,150))
lines(ss,res_randUneq[4,],typ="o",col=3,lwd=2,pch=15,ylim=c(0.5,1),xlim=c(40,150))
lines(ss,res_rand3Uneq[4,],typ="o",col=4,lwd=2,pch=15,ylim=c(0.5,1),xlim=c(40,150))
abline(h=c(0.8,0.9),lty=3,col=6,lwd=2)


plot(ss,res_sing[6,],typ="o",col=1,lwd=2,pch=15,ylim=c(0.3,1),xlim=c(40,150),main="alpha=0.2",
     ylab="Sample Siz",xlab="Power")
lines(ss,res_randEq[6,],typ="o",col=2,lwd=2,pch=15,ylim=c(0.5,1),xlim=c(40,150))
lines(ss,res_randUneq[6,],typ="o",col=3,lwd=2,pch=15,ylim=c(0.5,1),xlim=c(40,150))
lines(ss,res_rand3Uneq[6,],typ="o",col=4,lwd=2,pch=15,ylim=c(0.5,1),xlim=c(40,150))
abline(h=c(0.8,0.9),lty=3,col=6,lwd=2)




plot(ss,res_sing[4,],typ="0",col=1,lwd=2)

design_singArm2

### To tabulate
pwrEst(res_sing[1,],ss)
pwrEst(res_sing[2,],ss)
pwrEst(res_sing[3,],ss)
pwrEst(res_sing[4,],ss)
pwrEst(res_sing[5,],ss)
pwrEst(res_sing[6,],ss)

ss <- c(40,50,65,75,100,150)


### Figure to save
plot(ss,res[,1],typ="o",col=1,ylim=c(0.4,1),xlab="Sample Size",ylab="Power",pch=15)
points(ss,res[,2],typ="o",col=2,pch=15)
points(ss,res[,3],typ="o",col=3,pch=15)
points(ss,res[,4],typ="o",col=4,pch=15)
points(ss,res[,5],typ="o",col=5,pch=15)
points(ss,res[,6],typ="o",col=6,pch=15)
legend(130,0.75,c("0.01","0.025","0.05","0.1","0.15","0.2"),col=c(1:6),lty=1,bty="n")
abline(h=c(0.8,0.85,0.9),lty=3)




##############################################################################
##############################################################################
##############################################################################


##############################
##############################
############### Functions to be loaded first
##############################
##############################



### Function estimate power for each alpha level

pwrEst <- function(y,ss){

  md <- lm(ss~y+I(y^2))
  co <- coef(md)

  pw80 <- co[1]+co[2]*0.8+co[3]*0.8^2
  pw85 <- co[1]+co[2]*0.85+co[3]*0.85^2
  pw90 <- co[1]+co[2]*0.9+co[3]*0.9^2

  c(pw80,pw85,pw90)

}


CFM <- abCFM_prog

fuTime <- 12
rec <- rec1
n0=40
n1=40
nsim=4
nsim.psc=500
burn.psc=200
bound=0
direction="greater"
alpha_eval=c(0.01,0.025,0.05,0.1,0.15,0.2)

pscDesign <- function(CFM,n0=0,n1,beta,fuTime,recTime,rec=NULL,
                      nsim=4,nsim.psc=500,burn.psc=200,bound=0,
                      direction="greater",
                      alpha_eval=c(0.01,0.025,0.05,0.1,0.15,0.2)){

  #### Sorting out recruitment estimate
  if(is.null(rec)){
    n <- n0 + n1
    cumRec <- round(seq(0,n,length=recTime))
    mnthRec <- c(0,diff(cumRec))
    rec <- data.frame("SitesOpen"=1,"Monthly.Rec"=mnthRec,"Cumualtive.Rec."=cumRec)
  }

  if(!is.null(rec)){
    n <- max(rec$Cumualtive.Rec.)
  }

  ### Simulate Trials
  trialSim <- lapply(c(1:nsim),function(x){
    trialSamp(CFM=CFM,n0=n0,n1=n1,beta=beta,fuTime=fuTime,
              recTime=recTime,rec=rec,nsim.psc=nsim.psc,
              burn.psc=burn.psc)})


  ## Estimate proportion of posterior distribution > (<) bound
  trialEval <- lapply(trialSim,function(x){
    mn <- as.numeric(as.character(x$post_mn))
    sd <- as.numeric(as.character(x$post_sd))
    postEval(mn,sd,bound=bound,direction=direction)
  }
  )


  ## Evaluate bound against alpha level
  trialAlp <- lapply(trialEval,function(x) as.numeric(x<alpha_eval))
  pwrEst <- colSums(Reduce(rbind,trialAlp))/nsim

  ###
  estimate <- data.frame(alpha_eval,pwrEst)

  ### Sumamry of trial parameters
  trialSumm <- colMeans(Reduce(rbind,trialSim))

  ### Returning object
  ret <- list(trialSumm,estimate)
  return(ret)
}




### Function estimate power for each alpha level

pwrEst <- function(y,ss){

  md <- lm(ss~y+I(y^2))
  co <- coef(md)

  pw80 <- co[1]+co[2]*0.8+co[3]*0.8^2
  pw85 <- co[1]+co[2]*0.85+co[3]*0.85^2
  pw90 <- co[1]+co[2]*0.9+co[3]*0.9^2

  c(pw80,pw85,pw90)

}


##############################################################################
### Function for combinding estimates from psc and cm

pscComb <- function(x) {
  drs <- as_draws(x$draws)
  direct <- drs$beta_2
  indirect <- drs$beta_2 - drs$beta_1

  rho <- cor(direct, indirect)

  sig_d <- sd(direct)^2
  sig_i <- sd(indirect)^2

  mu_d <- mean(direct)
  mu_i <- mean(indirect)

  w <- sig_i / (sig_d + sig_i)

  mu_comb <- w * mu_d + (1 - w) * mu_i
  sig_comb <- w^2 *
    sig_d +
    (1 - w)^2 * sig_i +
    2 * w * (1 - w) * rho * sqrt(sig_d) * sqrt(sig_i)

  c(mu_comb, sqrt(sig_comb))
}


#######################################


pscfit <- function(
    CFM,
    DC,
    nsim = 2000,
    id = NULL,
    trt = NULL,
    nchain = 2,
    thin = 2,
    burn = 500
) {
  #### Step 1 - create pscCFM object (may not be required if pscCFM object supplied)
  if (!"pscCFM" %in% class(CFM)) {
    CFM <- pscCFM(CFM)
  }

  ### Step 2 - create pscfit data object
  pscOb <- pscData(CFM, DC, id = id, trt = trt)

  ### Step 3 - get initial values
  pscOb <- init(pscOb)

  ### Step 4 - MCMC estimation
  pscOb <- pscEst(pscOb, nsim, nchain)

  ### Step 5 - Summarising posterior
  pscOb <- postSummary(pscOb, thin = thin, burn = burn)

  ### Giving 'class to result'
  class(pscOb) <- "psc"
  pscOb
}


pscEst <- function(pscOb, nsim = 1000, nchain = 1) {
  ### Set Up
  pscOb <- pscEst_start(pscOb, nsim = nsim, nchain = nchain)
  pscOb

  ### Perform MCMC
  pscOb <- pscEst_run(pscOb, nsim = nsim, nchain = nchain)
  pscOb
}


pscEst_start <- function(pscOb, nsim, nchain) {
  #starting parameters
  if (!is.null(pscOb$trt)) {
    pscOb$start.mu <- rmvnorm(1, pscOb$start.mu, sigma = pscOb$start.sd)
  }

  if (is.null(pscOb$trt)) {
    pscOb$start.mu <- rnorm(1, pscOb$start.mu, pscOb$start.sd)
  }

  ### Setting up matrix to save draws
  draws <- matrix(NA, nsim, length(pscOb$co) + length(pscOb$start.mu) + 1)
  draws[1, ] <- c(pscOb$co, pscOb$start.mu, NA)
  beta.nm <- paste("beta", 1:length(pscOb$start.mu), sep = "_")
  colnames(draws) <- c(names(pscOb$co), beta.nm, "likEst")

  ## adding hazard parameters onto covariates for flexsurvsreg models
  cfmPost <- function(x) c(mvtnorm::rmvnorm(x, pscOb$co, pscOb$sig))

  #starting parameters
  if (!is.null(pscOb$trt)) {
    target <- function(x) {
      c(mvtnorm::rmvnorm(x, pscOb$start.mu, pscOb$start.sd * 2))
    }
  }

  if (is.null(pscOb$trt)) {
    target <- function(x) c(rnorm(x, pscOb$start.mu, pscOb$start.sd * 2))
  }

  # prior distributions
  betaPrior <- function(x) {
    mvtnorm::dmvnorm(x, rep(0, length(x)), diag(length(x)) * 1000, log = T)
  }

  ## ncores
  if (.Platform$OS.type == "windows" & nchain > 1) {
    warning("Currently only single chains allowed on Windows OS")
    nchain <- 1
  }

  # number of cores (for parallel computing of multiple chains)
  dcores <- detectCores()
  if (dcores > 3) {
    ncores <- min(round(dcores * .75), nchain)
  }


  ## Returning estimation object
  pscOb$cfmPost <- cfmPost
  pscOb$target <- target
  pscOb$betaPrior <- betaPrior

  pscOb$ncores <- ncores
  pscOb$draws <- draws


  ## Returning object
  pscOb
}




trialSamp <- function(CFM,n0,n1,beta,fuTime,recTime,rec,
                      nsim.psc=750,burn.psc=250){

  ##### Single arm estimation
  if(n0==0){
    ## Simualte data
    ds <- dataSim(CFM=CFM,n0=n0,n1=n1,beta=beta,recTime=recTime,
                  fuTime=fuTime,rec=rec)

    ## fitpsc
    simfit <- psc::pscfit(CFM,ds,nsim=nsim.psc,burn=burn.psc)

    co <- data.frame(coef(simfit))
    mn <- as.numeric(co$mean)
    sd <- as.numeric(co$sd)
    ret <- data.frame(sum(ds$cen),mn,sd)
    names(ret) <- c("ne","post_mn","post_sd")
  }



  #### Randomisation Estimation
  if(n0>0){

    ## Simualte data
    ds <- dataSim(CFM=CFM,n0=n0,n1=n1,beta=beta,recTime=recTime,
                  fuTime=fuTime,rec=rec)


    ### seperating dataset into arm 0 and arm 1
    #ds0 <- ds[ds$arm==0,]
    #ds1 <- ds[ds$arm==1,]

    ## fitpsc
    #simfit <- pscfit(CFM,ds1,nsim=nsim.psc)
    #co <- data.frame(coef(simfit))
    #mn_psc <- as.numeric(co$mean)
    #sd_psc <- as.numeric(co$sd);sd_psc

    ## direct estimation
    #cm <- coxph(Surv(ds$time,ds$cen)~ds$arm);cm
    #mn_d <- summary(cm)$coefficients[1]
    #sd_d <- summary(cm)$coefficients[3]

    ### if n0 is small cox model may not fit - and a warning is returned

    ### Bayesian Posterior Estimate
    #post_var <-  (1/sd_d^2+1/sd_psc^2)^(-1)
    #post_sd <- sqrt(post_var)
    #num <- ((mn_d/(sd_d^2))+(mn_psc/sd_psc^2))
    #post_mn <- num*post_var

    #ret <- data.frame(sum(ds$cen),mn_psc,sd_psc,mn_d,sd_d,post_mn,post_sd)
    #names(ret) <- c("ne","mn_psc","sd_psc","mn_direct","sd_direct","post_mn","post_sd")


    #### simfit2

    simfit <- pscfit(CFM,ds,nsim=nsim.psc,burn=burn.psc,trt=ds$arm,nchain=1)
    combRes <- pscComb(simfit)
    ret <- data.frame(sum(ds$cen),combRes[1],combRes[2])
    names(ret) <- c("ne","post_mn","post_sd")
    ret

  }

  ret
}





