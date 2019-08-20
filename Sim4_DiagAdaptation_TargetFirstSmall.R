## SET WORKING DIRECTORY AND LOAD FUNCTIONS
setwd('C:/Users/ZanellaG/Dropbox/BalancedMCMC/Balanced MCMC/Code_giacomo/Code for github/')
source('functions.R')

# define decaying schedule of learning rate
kappa<-0.6
gamma<-function(t,kappa_=kappa){
  return(t^(-kappa_))
}

# define target (100-dimensional Gaussian with first coordinate small)
n<-100 # number of dimensions
Sigma_targ<-diag(rep(1,n))
Sigma_targ[1,1]<-0.01^2
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

# DEFINE NUMBER OF MCMC ITERATIONS
T<-10^4

## RUN ALGORITHMS
output_MALA<-diag_grad_MCMC(T=T,n=n,method="MALA",start_x=start_x,start_sigma=start_sigma)
output_RWM<-diag_grad_MCMC(T=T,n=n,method="RWM",start_x=start_x,start_sigma=start_sigma)
output_Barker<-diag_grad_MCMC(T=T,n=n,method="Barker",start_x=start_x,start_sigma=start_sigma)

##############################
####### PLOTTING #############
##############################
par(mfrow=c(4,3))
time_window<-c(1:T)
ylim=range(c(output_MALA$sigma_vec,output_RWM$sigma_vec,output_Barker$sigma_vec))
plot(output_MALA$sigma_vec,type="l",ylim=ylim,main="MALA\n Adaptation of global scale",log = "y",ylab="",xlab="Markov chain iteration")
plot(output_RWM$sigma_vec,type="l",ylim=ylim,main="RWM\n Adaptation of global scale",log = "y",ylab="",xlab="Markov chain iteration")
plot(output_Barker$sigma_vec,type="l",ylim=ylim,main="Barker\n Adaptation of global scale",log = "y",ylab="",xlab="Markov chain iteration")
ylim<-NULL
for (i in c(1:n)){
  ylim=range(c(ylim,output_MALA$sigma_t[time_window,i],output_RWM$sigma_t[time_window,i],output_Barker$sigma_t[time_window,i]))
}
plot((output_MALA$sigma_t[time_window,3]),type="l",ylim=ylim,main="Adaptation of local scales",ylab="",xlab="Markov chain iteration",log="y")
for(j in c(1:n)){
  lines((output_MALA$sigma_t[time_window,j]))
}
plot((output_RWM$sigma_t[time_window,3]),type="l",ylim=ylim,main="Adaptation of local scales",ylab="",xlab="Markov chain iteration",log="y")
for(j in c(1:n)){
  lines((output_RWM$sigma_t[time_window,j]))
}
plot((output_Barker$sigma_t[time_window,3]),type="l",ylim=ylim,main="Adaptation of local scales",ylab="",xlab="Markov chain iteration",log="y")
for(j in c(1:n)){
  lines((output_Barker$sigma_t[time_window,j]))
}
i<-1
ylim=range(c(output_MALA$x_samples[time_window,i],output_RWM$x_samples[time_window,i],output_Barker$x_samples[time_window,i]))
ylim<-c(-0.1,0.1)
plot(output_MALA$x_samples[time_window,i],type="l",ylim=ylim,main="Traceplots of x_1",xlab="Markov chain iteration",ylab = "")
plot(output_RWM$x_samples[time_window,i],type="l",ylim=ylim,main="Traceplots of x_1",xlab="Markov chain iteration",ylab = "")
plot(output_Barker$x_samples[time_window,i],type="l",ylim=ylim,main="Traceplots of x_1",xlab="Markov chain iteration",ylab ="")
if(n>1){
  ylim<-NULL
  for (i in c(2:n)){
    ylim=range(c(ylim,output_MALA$x_samples[time_window,i],output_RWM$x_samples[time_window,i],output_Barker$x_samples[time_window,i]))
  }
  plot(output_MALA$x_samples[time_window,3],type="l",ylim=ylim,main="Traceplots of x_2,...,x_n",ylab ="",xlab="Markov chain iteration")
  for(j in c(2:n)){
    lines(output_MALA$x_samples[time_window,j])
  }
  plot(output_RWM$x_samples[time_window,3],type="l",ylim=ylim,main="Traceplots of x_2,...,x_n",ylab ="",xlab="Markov chain iteration")
  for(j in c(2:n)){
    lines(output_RWM$x_samples[time_window,j])
  }
  plot(output_Barker$x_samples[time_window,3],type="l",ylim=ylim,main="Traceplots of x_2,...,x_n",ylab ="",xlab="Markov chain iteration")
  for(j in c(2:n)){
    lines(output_Barker$x_samples[time_window,j])
  }
}
par(mfrow=c(1,1))