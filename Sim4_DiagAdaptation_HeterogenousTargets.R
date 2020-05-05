# the following code reproduces the simulations reported in 
# Section 6.2 "Performance on target distributions with heterogeneous scales"

## LOAD FUNCTIONS FROM GITHUB REPOSITORY
source('https://raw.githubusercontent.com/gzanella/barker/master/functions.R')

# define decaying schedule of learning rate
kappa<-0.6
gamma<-function(t,kappa_=kappa){
  return(t^(-kappa_))
}

# choose scenario under consideration (i.e. type of target distribution, see Section 6.2)
scenario<-1
# define corresponding target distribution (100-dimensional distribution with heterogneous scales)
set_target(scenario)

# CHOOSE RANDOM STARTING POINT
start_x<-rnorm(n,0,10)

# DEFINE NUMBER OF MCMC ITERATIONS
T<-10^4

## RUN ALGORITHMS
print("Running algorithms")
output_MALA<-diag_grad_MCMC(T=T,n=n,method="MALA",start_x=start_x,start_sigma=NULL)
output_RWM<-diag_grad_MCMC(T=T,n=n,method="RWM",start_x=start_x,start_sigma=NULL)
output_Barker<-diag_grad_MCMC(T=T,n=n,method="Barker",start_x=start_x,start_sigma=NULL)

## PLOT OUTPUT
print("Plotting output")
plot_tuning_traceplots(output_MALA,output_RWM,output_Barker)
