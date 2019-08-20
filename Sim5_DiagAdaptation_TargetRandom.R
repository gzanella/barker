## SET WORKING DIRECTORY ND LOAD PACKAGES
setwd('C:/Users/ZanellaG/Dropbox/BalancedMCMC/Balanced MCMC/Code_giacomo/Code for github/')
source('functions.R')

# define adaptation schedule
c1<-0.6
gamma_1<-function(t,c_0=1,c_1=c1){
  return(c_0*t^(-c_1))
}

###### DEFINE TARGET (100-dimensional Gaussian with random standard deviations)
n<-100
variances_sd<-3
variances<-exp(rnorm(n,mean = 0,sd = variances_sd))
Sigma_targ<-diag(variances)
Prec_target<-solve(Sigma_targ)
log_f_ratio<-function(x,y){return(-0.5*(
  matrix(y,nrow = 1,ncol = n)%*%Prec_target%*%matrix(y,nrow = n,ncol = 1)-
    matrix(x,nrow = 1,ncol = n)%*%Prec_target%*%matrix(x,nrow = n,ncol = 1)
))}
g_prime<-function(x){return(
  -c(matrix(x,nrow = 1,ncol = n)%*%Prec_target)
)}

# CHOOSE RANDOM STARTING POINT
start_x<-rnorm(n,0,10)
start_sigma<-NULL

# DEFINE NUMBER OF MONTE CARLO/MCMC ITERATIONS
T<-10^4

## RUN ALGORITHMS
output_MALA<-diag_grad_MCMC(T=T,n=n,method="MALA",start_x=start_x,start_sigma=start_sigma)
output_RWM<-diag_grad_MCMC(T=T,n=n,method="RWM",start_x=start_x,start_sigma=start_sigma)
output_Barker<-diag_grad_MCMC(T=T,n=n,method="Barker",start_x=start_x,start_sigma=start_sigma)

##############################
####### PLOTTING #############
##############################
par(mfrow=c(3,3))
time_window<-c(1:T)
ylim=range(c(output_MALA$sigma_vec,output_RWM$sigma_vec,output_Barker$sigma_vec))
plot(output_MALA$sigma_vec,type="l",ylim=ylim,main="MALA\n Adaptation of global scale",log = "y",ylab="",xlab="Markov chain iteration")
plot(output_RWM$sigma_vec,type="l",ylim=ylim,main="RWM\n Adaptation of global scale",log = "y",ylab="",xlab="Markov chain iteration")
plot(output_Barker$sigma_vec,type="l",ylim=ylim,main="Barker\n Adaptation of global scale",log = "y",ylab="",xlab="Markov chain iteration")
par(mar=c(2,2,1.5,0))
ylim<-NULL
for (i in c(1:n)){
  ylim=range(c(ylim,output_MALA$sigma_t[time_window,i]/sqrt(variances[i]),output_RWM$sigma_t[time_window,i]/sqrt(variances[i]),output_Barker$sigma_t[time_window,i]/sqrt(variances[i])))
}
plot((output_MALA$sigma_t[time_window,3]/sqrt(variances[3])),type="l",ylim=ylim,main="Adaptation of local scales (rescaled)",ylab="",xlab="Markov chain iteration",log="y")
for(j in c(1:n)){
  lines((output_MALA$sigma_t[time_window,j]/sqrt(variances[j])))
}
plot((output_RWM$sigma_t[time_window,3]/sqrt(variances[3])),type="l",ylim=ylim,main="Adaptation of local scales (rescaled)",ylab="",xlab="Markov chain iteration",log="y")
for(j in c(1:n)){
  lines((output_RWM$sigma_t[time_window,j]/sqrt(variances[j])))
}
plot((output_Barker$sigma_t[time_window,3]/sqrt(variances[3])),type="l",ylim=ylim,main="Adaptation of local scales (rescaled)",ylab="",xlab="Markov chain iteration",log="y")
for(j in c(1:n)){
  lines((output_Barker$sigma_t[time_window,j]/sqrt(variances[j])))
}
ylim<-NULL
for (i in c(1:n)){
  ylim=range(c(ylim,c(output_MALA$x_samples[time_window,i],output_RWM$x_samples[time_window,i],output_Barker$x_samples[time_window,i])/sqrt(variances[i])))
}
ylim[1]<-max(ylim[1],-60)
ylim[2]<-min(ylim[2],60)
plot(output_MALA$x_samples[time_window,3]/sqrt(variances[3]),type="l",ylim=ylim,main="Traceplots of x_1,...,x_n (rescaled)",ylab ="xi",xlab="Markov chain iteration")
for(j in c(1:n)){
  lines(output_MALA$x_samples[time_window,j]/sqrt(variances[j]))
}
plot(output_RWM$x_samples[time_window,3]/sqrt(variances[3]),type="l",ylim=ylim,main="Traceplots of x_1,...,x_n (rescaled)",ylab ="xi",xlab="Markov chain iteration")
for(j in c(1:n)){
  lines(output_RWM$x_samples[time_window,j]/sqrt(variances[j]))
}
plot(output_Barker$x_samples[time_window,3]/sqrt(variances[3]),type="l",ylim=ylim,main="Traceplots of x_1,...,x_n (rescaled)",ylab ="xi",xlab="Markov chain iteration")
for(j in c(1:n)){
  lines(output_Barker$x_samples[time_window,j]/sqrt(variances[j]))
}
par(mfrow=c(1,1),cex=1)