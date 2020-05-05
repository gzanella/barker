# the following code reproduces the simulation 
# reported in Section 5.1 "Illustrations of robustness to tuning"

## LOAD FUNCTIONS FROM GITHUB REPOSITORY
source('https://raw.githubusercontent.com/gzanella/barker/master/functions.R')

# DEFINE NUMBER OF MONTE CARLO/MCMC ITERATIONS
T<-10^3 # **** INCREASE T TO 10^5 TO GET MORE STABLE ESTIMATES ****

# DEFINE GRID OF VALUES OF SIGMA
sigma_vec<-2.38*exp(seq(from = -5,to = 6,length.out = 20))

### DEFINE TARGET (standard normal)
n<-1
log_f_ratio<-function(x,y){return(sum(-y^2/2+x^2/2))}
g_prime<-function(x){return(-x)}
rtarget<-function(n){return(rnorm(1))}

### COMPUTE ESJD OVER STEP-SIZE FOR THE THREE SCHEMES
ESJD_MALA_G<-rep(NA,length(sigma_vec))
for (sigma in sigma_vec){
  ESJD_MALA_G[which(sigma_vec==sigma)]<-grad_MCMC_ESJD(T=T,n=n,method = "MALA",sigma=sigma)
}
ESJD_RW_G<-rep(NA,length(sigma_vec))
for (sigma in sigma_vec){
  ESJD_RW_G[which(sigma_vec==sigma)]<-grad_MCMC_ESJD(T=T,n=n,method = "RWM",sigma=sigma)
}
ESJD_BARKER_G<-rep(NA,length(sigma_vec))
for (sigma in sigma_vec){
  ESJD_BARKER_G[which(sigma_vec==sigma)]<-grad_MCMC_ESJD(T=T,n=n,method = "Barker",sigma=sigma)
}

### DEFINE TARGET (laplace) AND GRID OF VALUES OF SIGMA ####
log_f_ratio<-function(x,y){return( -abs(y)+abs(x) ) }
g_prime<-function(x){return(-sign(x) )}
rtarget<-function(n){return(sign(runif(n)-0.5)*rexp(n=n))}
### COMPUTE ESJD OVER STEP-SIZE FOR THE THREE SCHEMES
ESJD_MALA_L<-rep(NA,length(sigma_vec))
for (sigma in sigma_vec){
  ESJD_MALA_L[which(sigma_vec==sigma)]<-grad_MCMC_ESJD(T=T,n=n,method = "MALA",sigma=sigma)
}
ESJD_RW_L<-rep(NA,length(sigma_vec))
for (sigma in sigma_vec){
  ESJD_RW_L[which(sigma_vec==sigma)]<-grad_MCMC_ESJD(T=T,n=n,method = "RWM",sigma=sigma)
}
ESJD_BARKER_L<-rep(NA,length(sigma_vec))
for (sigma in sigma_vec){
  ESJD_BARKER_L[which(sigma_vec==sigma)]<-grad_MCMC_ESJD(T=T,n=n,method = "Barker",sigma=sigma)
}

### DEFINE TARGET (student-t) AND GRID OF VALUES OF SIGMA ####
df<-5
log_f_ratio<-function(x,y){return(-(df+1)/2*log(1+y^2/df)+(df+1)/2*log(1+x^2/df))}
g_prime<-function(x){return(-(df+1)*x/(df+x^2) )}
rtarget<-function(n){return(rt(n=n,df=df))}
### COMPUTE ESJD OVER STEP-SIZE FOR THE THREE SCHEMES
ESJD_MALA_T<-rep(NA,length(sigma_vec))
for (sigma in sigma_vec){
  ESJD_MALA_T[which(sigma_vec==sigma)]<-grad_MCMC_ESJD(T=T,n=n,method = "MALA",sigma=sigma)
}
ESJD_RW_T<-rep(NA,length(sigma_vec))
for (sigma in sigma_vec){
  ESJD_RW_T[which(sigma_vec==sigma)]<-grad_MCMC_ESJD(T=T,n=n,method = "RWM",sigma=sigma)
}
ESJD_BARKER_T<-rep(NA,length(sigma_vec))
for (sigma in sigma_vec){
  ESJD_BARKER_T[which(sigma_vec==sigma)]<-grad_MCMC_ESJD(T=T,n=n,method = "Barker",sigma=sigma)
}

##############################
####### PLOTTING #############
##############################

# PLOT ESJD OF FIRST COORDINATE
par(mfrow=c(1,3))
par(mar=c(3,2.7,1.6,1))
# GAUSSIAN TARGET
log_axes<-"xy"
ylim=c(exp(-13),max(ESJD_MALA_G,ESJD_RW_G,ESJD_BARKER_G,na.rm = TRUE))
plot(sigma_vec,ESJD_MALA_G,log=log_axes,col="red",
     xlab="proposal step-size",ylab="ESJD",ylim=ylim,pch=17)
lines(sigma_vec,ESJD_MALA_G,col="red")
points(sigma_vec,ESJD_RW_G,col="black",pch=16)
lines(sigma_vec,ESJD_RW_G,col="black")
points(sigma_vec,ESJD_BARKER_G,col="blue",pch=18)
lines(sigma_vec,ESJD_BARKER_G,col="blue")
title("Gaussian target")
legend(x = "bottomleft",legend = c("RW","MALA","Barker"),col = c("black","red","blue"),lwd = c(1,1,1),pch = c(16,17,18))
# LAPLACE TARGET
ylim=c(exp(-13),max(ESJD_MALA_L,ESJD_RW_L,ESJD_BARKER_L))
plot(sigma_vec,ESJD_MALA_L,log=log_axes,col="red",
     xlab="proposal step-size",ylab="ESJD",ylim=ylim,pch=17)
lines(sigma_vec,ESJD_MALA_L,col="red")
points(sigma_vec,ESJD_RW_L,col="black",pch=16)
lines(sigma_vec,ESJD_RW_L,col="black")
points(sigma_vec,ESJD_BARKER_L,col="blue",pch=18)
lines(sigma_vec,ESJD_BARKER_L,col="blue")
title("Laplace target")
# STUDENT-T TARGET
par(mar=c(3,2.7,1.6,1))
ylim=c(exp(-13),max(ESJD_MALA_T,ESJD_RW_T,ESJD_BARKER_T))
plot(sigma_vec,ESJD_MALA_T,log=log_axes,col="red",
     xlab="proposal step-size",ylab="ESJD",ylim=ylim,pch=17)
lines(sigma_vec,ESJD_MALA_T,col="red")
points(sigma_vec,ESJD_RW_T,col="black",pch=16)
lines(sigma_vec,ESJD_RW_T,col="black")
points(sigma_vec,ESJD_BARKER_T,col="blue",pch=18)
lines(sigma_vec,ESJD_BARKER_T,col="blue")
title(paste("Student-t target (df=",df,")",sep=""))
par(mfrow=c(1,1))
