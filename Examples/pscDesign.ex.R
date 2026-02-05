

### HCC Example

#install.packages("pscDesign")

library(psc)
library(pscDesign)


### Getting Atezo Bev model

setwd("/Users/richardjackson/Documents/GitHub/synthetic-Controlled-Trials/Examples/Data")

load("abCFM_prog.Rda")


### simulation paramteres (Once all scenarios complete can you increase these to 500, 2000, 400)
nsim = 50,
nsim.psc = 500,
burn.psc = 200,


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


#######################################
### single Arm design
design_singArm1 <- pscDesign(abCFM_prog,beta=beta,n0=0,n1=40,fuTime=12,rec=rec1,nsim=nsim)
design_singArm2 <- pscDesign(abCFM_prog,beta=beta,n0=0,n1=50,fuTime=12,rec=rec2,nsim=nsim)
design_singArm3 <- pscDesign(abCFM_prog,beta=beta,n0=0,n1=65,fuTime=12,rec=rec3,nsim=nsim)
design_singArm4 <- pscDesign(abCFM_prog,beta=beta,n0=0,n1=75,fuTime=12,rec=rec4,nsim=nsim)
design_singArm5 <- pscDesign(abCFM_prog,beta=beta,n0=0,n1=100,fuTime=12,rec=rec5,nsim=nsim)
design_singArm6 <- pscDesign(abCFM_prog,beta=beta,n0=0,n1=150,fuTime=12,rec=rec6,nsim=nsim)


#############################. Summarising results

res <- cbind(design_singArm1[[2]][,2],design_singArm2[[2]][,2],design_singArm3[[2]][,2],
             design_singArm4[[2]][,2],design_singArm5[[2]][,2],design_singArm6[[2]][,2])


### To tabulate
pwrEst(res[,1],ss)
pwrEst(res[,2],ss)
pwrEst(res[,3],ss)
pwrEst(res[,4],ss)
pwrEst(res[,5],ss)
pwrEst(res[,6],ss)

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


#############################


### Randomised (1:1)
design_randEq1 <- pscDesign(abCFM_prog,beta=beta,n0=20,n1=20,fuTime=12,rec=rec1,nsim=nsim)
design_randEq2 <- pscDesign(abCFM_prog,beta=beta,n0=25,n1=25,fuTime=12,rec=rec2,nsim=nsim)
design_randEq3 <- pscDesign(abCFM_prog,beta=beta,n0=32,n1=33,fuTime=12,rec=rec3,nsim=nsim)
design_randEq4 <- pscDesign(abCFM_prog,beta=beta,n0=37,n1=38,fuTime=12,rec=rec4,nsim=nsim)
design_randEq5 <- pscDesign(abCFM_prog,beta=beta,n0=50,n1=50,fuTime=12,rec=rec5,nsim=nsim)
design_randEq6 <- pscDesign(abCFM_prog,beta=beta,n0=075,n1=150,fuTime=12,rec=rec6,nsim=nsim)

### Randomised (1:2)
design_randUneq1 <- pscDesign(abCFM_prog,beta=beta,n0=20,n1=20,fuTime=12,rec=rec1,nsim=nsim)
design_randUneq2 <- pscDesign(abCFM_prog,beta=beta,n0=25,n1=25,fuTime=12,rec=rec2,nsim=nsim)
design_randUneq3 <- pscDesign(abCFM_prog,beta=beta,n0=32,n1=33,fuTime=12,rec=rec3,nsim=nsim)
design_randUneq4 <- pscDesign(abCFM_prog,beta=beta,n0=37,n1=38,fuTime=12,rec=rec4,nsim=nsim)
design_randUneq5 <- pscDesign(abCFM_prog,beta=beta,n0=50,n1=50,fuTime=12,rec=rec5,nsim=nsim)
design_randUneq6 <- pscDesign(abCFM_prog,beta=beta,n0=075,n1=150,fuTime=12,rec=rec6,nsim=nsim)



### Repeat with inflated variance
abCFM_prog_iv <- abCFM_prog
abm_sig <- abCFM_prog_iv$sig
diag(abm_sig) <- diag(abm_sig)*2
abCFM_prog_iv$sig <- abm_sig

### Repeat the 3 designs above with inflated variance

### Repeat all 6 designs with gemCFM model


###### For both gemCFM and abCFM_prog: For a given sample size of 65 (rec3)
###### What allocation ratio gives the most power?






### Function estimate power for each alpha level

pwrEst <- function(y,ss){

  md <- lm(ss~y+I(y^2))
  co <- coef(md)

  pw80 <- co[1]+co[2]*0.8+co[3]*0.8^2
  pw85 <- co[1]+co[2]*0.85+co[3]*0.85^2
  pw90 <- co[1]+co[2]*0.9+co[3]*0.9^2

  c(pw80,pw85,pw90)

}

design_singArm1

pscDesign

design_singArm

