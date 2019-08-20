setwd('C:/Users/ZanellaG/Dropbox/BalancedMCMC/Balanced MCMC/Code_giacomo/Code for github/')
source('functions.R')

### DEFINE ISOTROPIC NORMAL TARGET ####
log_f_ratio<-function(x,y){return(sum(-y^2/2+x^2/2))}
g_prime<-function(x){return(-x)}
rtarget<-function(n){return(rnorm(n))}
### T IS THE NUMBER OF MONTE CARLO SAMPLES ####
T<-1000

### n_vec IS THE GRID OF DIMENSIONS TO PLOT ####
n_vec<-ceiling(2^c(3:10))

## COMPUTE ESJD FOR EACH DIMENSION AND THEN TAKE THE MAXIMUM OVER SIGMA
max_ESJD_MALA<-rep(NA,length(n_vec))
max_ESJD_RW<-rep(NA,length(n_vec))
max_ESJD_BARKER<-rep(NA,length(n_vec))
for (n in n_vec){
  print(paste("n=",n))
  # choose grid of values for sigma around the theoretically optimal ones (NB: if needed use additional plots to check that grid is appropriate)
  sigma_vec<-2.38/sqrt(n)*exp(seq(from = log(0.4),to = log(3),length.out = 20))
  sigma_vec_grad<-2.38/(n^(1/6))*exp(seq(from = log(0.1),to = log(1.5),length.out = 20))
  ## COMPUTE ESJD OVER STEP-SIZE FOR THE THREE SCHEMES AND TAKE THE OPTIMAL ONE FOR EACH SCHEME
  ### MALA ####
  ESJD_MALA<-rep(NA,length(sigma_vec_grad))
  for (sigma in sigma_vec_grad){
    ESJD_MALA[which(sigma_vec_grad==sigma)]<-mean(grad_MCMC_ESJD(T=T,n=n,method = "MALA",sigma = sigma))
  }
  ### RWM ####
  ESJD_RW<-rep(NA,length(sigma_vec))
  for (sigma in sigma_vec){
    ESJD_RW[which(sigma_vec==sigma)]<-mean(grad_MCMC_ESJD(T=T,n=n,method = "RWM",sigma = sigma))
  }
  ### BARKER ####
  ESJD_BARKER<-rep(NA,length(sigma_vec_grad))
  for (sigma in sigma_vec_grad){
    ESJD_BARKER[which(sigma_vec_grad==sigma)]<-mean(grad_MCMC_ESJD(T=T,n=n,method = "Barker",sigma = sigma))
  }
  iter<-which(n==n_vec)
  max_ESJD_MALA[iter]<-max(ESJD_MALA)
  max_ESJD_RW[iter]<-max(ESJD_RW)
  max_ESJD_BARKER[iter]<-max(ESJD_BARKER)
}


##############################
####### PLOTTING #############
##############################
## PLOT ESJD OVER DIMENSION
log_axis<-"xy"
plot(n_vec,max_ESJD_MALA,col="red",log=log_axis,ylab="ESJD per coordinate",xlab="number of dimensions",
     ylim=range(max_ESJD_MALA,max_ESJD_RW,max_ESJD_BARKER),
     pch=17)
lines(n_vec,max_ESJD_MALA,col="red")
lines((n_vec),max_ESJD_RW)
points((n_vec),max_ESJD_RW,pch=16)
lines((n_vec),max_ESJD_BARKER,col="blue")
points((n_vec),max_ESJD_BARKER,col="blue",pch=18)
legend(x = "bottomleft",legend = c("RW","MALA","Barker"),col = c("black","red","blue"),lwd = c(1,1,1),pch = c(16,17,18))
title("Efficiency over dimensionality - isotropic target")