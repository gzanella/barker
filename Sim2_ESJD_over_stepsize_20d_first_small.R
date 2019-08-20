setwd('C:/Users/ZanellaG/Dropbox/BalancedMCMC/Balanced MCMC/Code_giacomo/Code for github/')
source('functions.R')

# ### DEFINE TARGET (MULTIV-NORMAL WITH FIRST COMPONENT SMALL) ####
n<-20
sigma_targets<-c(0.01,rep(1,n-1))
log_f_ratio<-function(x,y,sigma_targ=sigma_targets){return(sum(-(y/sigma_targ)^2/2+(x/sigma_targ)^2/2))}
g_prime<-function(x,sigma_targ=sigma_targets){return(-x/(sigma_targ^2))}
rtarget<-function(n,sigma_targ=sigma_targets){return(rnorm(n,mean = 0,sd = sigma_targ))}

# DEFINE NUMBER OF MCMC ITERATIONS
T<-10^4 # !!! INCREASE T TO 10^5 TO GET MORE STABLE ESTIMATES !!!

# DEFINE GRID OF VALUES OF SIGMA
sigma_vec<-2.38/sqrt(n)*exp(seq(from = -7.5,to = 2,length.out = 25))

## COMPUTE ESJD OVER STEP-SIZE FOR THE THREE SCHEMES
### MALA ####
ESJD_MALA<-matrix(rep(NA,length(sigma_vec)*n),nrow = length(sigma_vec),ncol = n)
for (sigma in sigma_vec){
  ESJD_MALA[which(sigma_vec==sigma),]<-grad_MCMC_ESJD(T=T,n=n,method = "MALA",sigma=sigma)
}
ESJD_MALA<-ESJD_MALA#/rep(sigma_targets,each=length(sigma_vec))# normalize for target sd
### RWM ####
ESJD_RW<-matrix(rep(NA,length(sigma_vec)*n),nrow = length(sigma_vec),ncol = n)
for (sigma in sigma_vec){
  ESJD_RW[which(sigma_vec==sigma),]<-grad_MCMC_ESJD(T=T,n=n,method = "RWM",sigma=sigma)
}
ESJD_RW<-ESJD_RW#/rep(sigma_targets,each=length(sigma_vec))# normalize for target sd
### BARKER ####
ESJD_BARKER<-matrix(rep(NA,length(sigma_vec)*n),nrow = length(sigma_vec),ncol = n)
for (sigma in sigma_vec){
  ESJD_BARKER[which(sigma_vec==sigma),]<-grad_MCMC_ESJD(T=T,n=n,method = "Barker",sigma=sigma)
}
ESJD_BARKER<-ESJD_BARKER#/rep(sigma_targets,each=length(sigma_vec))# normalize for target sd

##############################
####### PLOTTING #############
##############################
# PLOT ESJD OF FIRST COORDINATE AND MEDIAN ESJD OF OTHER COORDINATES
par(mfrow=c(1,2))
#plot ESJD of first coordinate
log_axes<-"xy"
ylim=c(exp(-15),max(ESJD_MALA[,1],ESJD_RW[,1],ESJD_BARKER[,1]))
plot(sigma_vec,ESJD_MALA[,1],log=log_axes,col="red",
     xlab="proposal step-size",ylab="ESJD",ylim=ylim,pch=17)
lines(sigma_vec,ESJD_MALA[,1],col="red")
points(sigma_vec,ESJD_RW[,1],col="black",pch=16)
lines(sigma_vec,ESJD_RW[,1],col="black")
points(sigma_vec,ESJD_BARKER[,1],col="blue",pch=18)
lines(sigma_vec,ESJD_BARKER[,1],col="blue")
title("ESJD of coordinate 1")
legend(x = "topleft",legend = c("RW","MALA","Barker"),col = c("black","red","blue"),lwd = c(1,1,1),pch = c(16,17,18))
#plot median ESJD OF OTHER COORDINATES
log_axes<-"xy"
mala_median<-apply(ESJD_MALA[,-1],MARGIN = 1,median)
rwm_median<-apply(ESJD_RW[,-1],MARGIN = 1,median)
barker_median<-apply(ESJD_BARKER[,-1],MARGIN = 1,median)
ylim=c(exp(-15),max(mala_median,rwm_median,barker_median))
plot(sigma_vec,mala_median,log=log_axes,col="red",
     xlab="proposal step-size",ylab="ESJD",ylim=ylim,pch=17)
lines(sigma_vec,mala_median,col="red")
points(sigma_vec,rwm_median,col="black",pch=16)
lines(sigma_vec,rwm_median,col="black")
points(sigma_vec,barker_median,col="blue",pch=18)
lines(sigma_vec,barker_median,col="blue")
title("median ESJD of coordinates 2 to 20")
par(mfrow=c(1,1))