

### HCC Example
rm(list=ls())
setwd("/Users/richardjackson/Documents/GitHub/psc")
devtools::load_all()
setwd("/Users/richardjackson/Documents/GitHub/pscDesign")
devtools::load_all()

#library(psc)
#library(pscDesign)


### Getting Atezo Bev model

setwd("/Users/richardjackson/Documents/GitHub/synthetic-Controlled-Trials/Examples/Data")
dir()

load("abCFM_prog_ex.Rda")
load("abCFM_prog_ex2.Rda")
load("abCFM_prog.Rda")





### simulation parameters (Once all scenarios complete can you increase these to 500, 2000, 400)
nsim = 1000
nsim.psc = 5000
burn.psc = 500

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

recForcast

rec.forcast <- recForcast


#######################################
setwd("SimRes")
#######################################
### single Arm design
design_singArm1 <- pscDesign(abCFM_prog,beta=log(0.725),n0=0,n1=40,fuTime=12,rec=rec1,nsim=nsim,cen.prop=0.5)
save(design_singArm1,file="design_singArm1.2.R")
design_singArm2 <- pscDesign(abCFM_prog,beta=log(0.725),n0=0,n1=50,fuTime=12,rec=rec2,nsim=nsim,cen.prop=0.5)
save(design_singArm2,file="design_singArm2.2.R")
design_singArm3 <- pscDesign(abCFM_prog,beta=log(0.725),n0=0,n1=65,fuTime=12,rec=rec3,nsim=nsim,cen.prop=0.5)
save(design_singArm3,file="design_singArm3.2.R")
design_singArm4 <- pscDesign(abCFM_prog,beta=log(0.725),n0=0,n1=75,fuTime=12,rec=rec4,nsim=nsim,cen.prop=0.5)
save(design_singArm4,file="design_singArm4.2.R")
design_singArm5 <- pscDesign(abCFM_prog,beta=log(0.725),n0=0,n1=100,fuTime=12,rec=rec5,nsim=nsim,cen.prop=0.5)
save(design_singArm5,file="design_singArm5.2.R")
design_singArm6 <- pscDesign(abCFM_prog,beta=log(0.725),n0=0,n1=150,fuTime=12,rec=rec6,nsim=nsim,cen.prop=0.5)
save(design_singArm6,file="design_singArm6.2.R")

#######################################
### Randomised (1:1)
design_randEq1 <- pscDesign(abCFM_prog,beta=log(0.725),n0=20,n1=20,fuTime=12,rec=rec1,nsim=nsim,cen.prop=0.5)
save(design_randEq1,file="design_randEq1.2.R")
design_randEq2 <- pscDesign(abCFM_prog,beta=log(0.725),n0=25,n1=25,fuTime=12,rec=rec2,nsim=nsim,cen.prop=0.5)
save(design_randEq2,file="design_randEq2.2.R")
design_randEq3 <- pscDesign(abCFM_prog,beta=log(0.725),n0=32,n1=33,fuTime=12,rec=rec3,nsim=nsim,cen.prop=0.5)
save(design_randEq3,file="design_randEq3.2.R")
design_randEq4 <- pscDesign(abCFM_prog,beta=log(0.725),n0=37,n1=38,fuTime=12,rec=rec4,nsim=nsim,cen.prop=0.5)
save(design_randEq4,file="design_randEq4.2.R")
design_randEq5 <- pscDesign(abCFM_prog,beta=log(0.725),n0=50,n1=50,fuTime=12,rec=rec5,nsim=nsim,cen.prop=0.5)
save(design_randEq5,file="design_randEq5.2.R")
design_randEq6 <- pscDesign(abCFM_prog,beta=log(0.725),n0=75,n1=75,fuTime=12,rec=rec6,nsim=nsim,cen.prop=0.5)
save(design_randEq6,file="design_randEq6.2.R")


### Randomised (1:2)
design_randUneq1 <- pscDesign(abCFM_prog,beta=log(0.725),n0=13,n1=27,fuTime=12,rec=rec1,nsim=nsim,cen.prop=0.5)
save(design_randUneq1,file="design_randUneq1.2.R")
design_randUneq2 <- pscDesign(abCFM_prog,beta=log(0.725),n0=17,n1=33,fuTime=12,rec=rec2,nsim=nsim,cen.prop=0.5)
save(design_randUneq2,file="design_randUneq2.2.R")
design_randUneq3 <- pscDesign(abCFM_prog,beta=log(0.725),n0=22,n1=43,fuTime=12,rec=rec3,nsim=nsim,cen.prop=0.5)
save(design_randUneq3,file="design_randUneq3.2.R")
design_randUneq4 <- pscDesign(abCFM_prog,beta=log(0.725),n0=25,n1=50,fuTime=12,rec=rec4,nsim=nsim,cen.prop=0.5)
save(design_randUneq4,file="design_randUneq4.2.R")
design_randUneq5 <- pscDesign(abCFM_prog,beta=log(0.725),n0=33,n1=67,fuTime=12,rec=rec5,nsim=nsim,cen.prop=0.5)
save(design_randUneq5,file="design_randUneq5.2.R")
design_randUneq6 <- pscDesign(abCFM_prog,beta=log(0.725),n0=50,n1=100,fuTime=12,rec=rec6,nsim=nsim,cen.prop=0.5)
save(design_randUneq6,file="design_randUneq6.2.R")

### Randomised (1:3)
design_rand3Uneq1 <- pscDesign(abCFM_prog,beta=log(0.725),n0=10,n1=30,fuTime=12,rec=rec1,nsim=nsim,cen.prop=0.5)
save(design_rand3Uneq1,file="design_rand3Uneq1.2.R")
design_rand3Uneq2 <- pscDesign(abCFM_prog,beta=log(0.725),n0=12,n1=38,fuTime=12,rec=rec2,nsim=nsim,cen.prop=0.5)
save(design_rand3Uneq2,file="design_rand3Uneq2.2.R")
design_rand3Uneq3 <- pscDesign(abCFM_prog,beta=log(0.725),n0=16,n1=49,fuTime=12,rec=rec3,nsim=nsim,cen.prop=0.5)
save(design_rand3Uneq3,file="design_rand3Uneq3.2.R")
design_rand3Uneq4 <- pscDesign(abCFM_prog,beta=log(0.725),n0=19,n1=56,fuTime=12,rec=rec4,nsim=nsim,cen.prop=0.5)
save(design_rand3Uneq4,file="design_rand3Uneq4.2.R")
design_rand3Uneq5 <- pscDesign(abCFM_prog,beta=log(0.725),n0=25,n1=75,fuTime=12,rec=rec5,nsim=nsim,cen.prop=0.5)
save(design_rand3Uneq5,file="design_rand3Uneq5.2.R")
design_rand3Uneq6 <- pscDesign(abCFM_prog,beta=log(0.725),n0=30,n1=120,fuTime=12,rec=rec6,nsim=nsim,cen.prop=0.5)
save(design_rand3Uneq6,file="design_rand3Uneq6.2.R")

getwd()

#############################. Summarising results. (TABLE)


design_singArm5


setwd("/Users/richardjackson/Documents/GitHub/synthetic-Controlled-Trials/Examples/Data/SimRes")

pwEst <- function(ds,a){
  tb <- table(ds$'P>0'>a)
  ptb <- tb/sum(tb)
  ptb[1]
}


### reading in
for(i in dir()) load(i)

alp <- c(0.05,0.1,0.2)

sa_des <- list(design_singArm1,design_singArm2,design_singArm3,
               design_singArm4,design_singArm5,design_singArm6)
ra_des <- list(design_randEq1,design_randEq2,design_randEq3,
               design_randEq4,design_randEq5,design_randEq6)
ura_des <- list(design_randUneq1,design_randUneq2,design_randUneq3,
               design_randUneq4,design_randUneq5,design_randUneq6)
ura3_des <- list(design_rand3Uneq1,design_rand3Uneq2,design_rand3Uneq3,
               design_rand3Uneq4,design_rand3Uneq5,design_rand3Uneq6)

sa_pow <- lapply(sa_des,function(y) lapply(alp,function(x) pwEst(y,x)));sa_pow
ra_pow <- lapply(ra_des,function(y) lapply(alp,function(x) pwEst(y,x)))
ura_pow <- lapply(ura_des,function(y) lapply(alp,function(x) pwEst(y,x)))
ura3_pow <- lapply(ura3_des,function(y) lapply(alp,function(x) pwEst(y,x)))

pw <- rep(c(0.05,0.1,0.2),6)
ss <- c(40,50,65,75,100,150)[gl(6,3)]


sa_pow[[5]]
ra_pow[[5]]
ura_pow[[5]]
ura3_pow[[5]]


exp(mean(design_rand3Uneq5$mean))
mean(design_rand3Uneq5$mean)
mean(design_rand3Uneq5$sd)

mean(design_randUneq5$mean)
mean(design_randUneq5$sd)

mean(design_randEq5$mean)
mean(design_randEq5$sd)



pwr_res <- data.frame(pw,ss,"sa"=unlist(sa_pow),"ra"=unlist(ra_pow),"ura"=unlist(ura_pow),"ura3"=unlist(ura3_pow))

library(tidyr)
pwr_lng <- pivot_longer(pwr_res,cols=3:6)

factor(pwr_lng$name)

ggplot(pwr_lng,aes(x=ss,y=value,group=name,color=name))+
  geom_line(linewidth=1.3)+
  facet_wrap(~pw)


####
k <- c(3,4,6)
pow <- c(0.8,0.85,0.9)

res_1 <- NULL
res_2 <- NULL
res_3 <- NULL
res_4 <- NULL
ss <- c(40,50,65,75,100,150)
res_sing

m<-3
n<-1
for(n in 1:3){
  for(m in 1:3){
    ln <- k[m]
    pw <- pow[n]
    ap <- spline(res_sing[ln,],ss,xout=pw);ap
    ap2 <- spline(res_sing[ln,],res_sing[7,],xout=pw);ap2
    re <- paste(round(ap$y)," (",round(ap2$y),")",sep="")
    res_1 <- c(res_1,re)

    ap <- spline(res_randEq[ln,],ss,xout=pw);ap
    ap2 <- spline(res_randEq[ln,],res_randEq[7,],xout=pw);ap2
    re <- paste(round(ap$y)," (",round(ap2$y),")",sep="")
    res_2 <- c(res_2,re)

    ap <- spline(res_randUneq[ln,],ss,xout=pw);ap
    ap2 <- spline(res_randUneq[ln,],res_randUneq[7,],xout=pw);ap2
    re <- paste(round(ap$y)," (",round(ap2$y),")",sep="")
    res_3 <- c(res_3,re)

    ap <- spline(res_rand3Uneq[ln,],ss,xout=pw);ap
    ap2 <- spline(res_rand3Uneq[ln,],res_rand3Uneq[7,],xout=pw);ap2
    re <- paste(round(ap$y)," (",round(ap2$y),")",sep="")
    res_4 <- c(res_4,re)
  }
}


sd





##############################################################################
##############################################################################
##############################################################################


##############################
##############################
############### Functions to be loaded first
##############################
##############################

CFM <- abCFM_prog
beta=log(0.725)
n0=0
n1=40
fuTime=12
rec=rec1
nsim=2
cen.prop=0.5
nsim=4
nsim.psc=500
burn.psc=200
bound=0
direction="greater"
alpha_eval=c(0.01,0.025,0.05,0.1,0.15,0.2)

pscDesign <- function(CFM,n0=0,n1,beta,fuTime,recTime,rec=NULL,
                      nsim=4,nsim.psc=500,burn.psc=200,bound=0,
                      direction="greater",
                      alpha_eval=c(0.01,0.025,0.05,0.1,0.15,0.2),
                      cen.prop=0.05){

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
  trialSim <- lapply(1:nsim,function(x){
    trialSamp(CFM=CFM,beta=beta,n0=n0,n1=n1,fuTime=fuTime,rec=rec,
              nsim.psc=nsim.psc,burn.psc=burn.psc,cen.prop=cen.prop)})

  trialSamp(CFM=CFM,beta=beta,n0=n0,n1=n1,fuTime=fuTime,rec=rec,
            nsim.psc=nsim.psc,burn.psc=burn.psc,cen.prop=cen.prop)
  trialSamp

  ### Organise Results
  l1 <- unlist(lapply(trialSim,mean))
  l2 <- unlist(lapply(trialSim,sd))
  l3 <- lapply(trialSim,quantile,p=c(0.025,0.5,0.975))
  l4 <- unlist(lapply(trialSim,postEval,bound=bound))

  ### Return
  res <- data.frame(cbind(l1,l2,Reduce("rbind",l3),l4))
  names(res)  <- c("mean","sd","2.5%","50%","97.5%","P>0")

  res

}





