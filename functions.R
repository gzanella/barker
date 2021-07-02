# this file contains the R functions used by the other main files to 
# reproduce the simulations reported in the paper

# required packages
require(shapes)
require(mvtnorm)

###################################################
####### ISOTROPIC PROPOSALS - NO ADAPTATION #######
###################################################
# this function implements RWM/MALA/Barker with isotropic proposal and no adaptation
grad_MCMC<-function(T=T,sigma=sigma,n=n,method=NULL,start_x){
  stopifnot(method=="MALA"|method=="RWM"|method=="Barker")
  # specify proposal and acceptance probability
  if(method=="MALA"){
    log_q_ratio<-log_q_ratio_mala
    rprop<-rmala
  }
  if(method=="RWM"){
    log_q_ratio<-log_q_ratio_rw
    rprop<-rrw
  }
  if(method=="Barker"){
    log_q_ratio<-log_q_ratio_barker
    rprop<-rbarker
  }
  # initialize sampler and output objects
  x<-start_x
  t<-1
  x_samples<-matrix(NA,nrow = T,ncol = n)
  x_samples[1,]<-x
  # mcmc iterations
  for (t in 2:T){
    # propose new state
    y<-x+rprop(g_prime(x),sigma)
    # compute acceptance rate
    ap<-min(1,exp(log_f_ratio(x,y)+log_q_ratio(x,y,sigma)))
    if (runif(1)<ap){# accept/reject
      x<-y
    }
    # store samples
    x_samples[t,]<-x
  }
  return(list(x_samples=x_samples))
}
## MALA PROPOSAL AND ACCEPTANCE PROBABILITY (ISOTROPIC)
log_q_ratio_mala<-function(x,y,sigma){return(sum(
  -(x-y-sigma^2*g_prime(y)/2)^2/(2*sigma^2)+
    (y-x-sigma^2*g_prime(x)/2)^2/(2*sigma^2)
))}
rmala <-function(c.,sigma.){
  return(
    sigma.*rnorm(n=length(c.))+sigma.^2*c./2
  )}
## RANDOM WALK PROPOSAL AND ACCEPTANCE PROBABILITY (ISOTROPIC)
log_q_ratio_rw<-function(x,y,sigma){return(0)}
rrw<-function(c,sigma){return(sigma*rnorm(n=length(c)))}
## BARKER PROPOSAL AND ACCEPTANCE PROBABILITY (ISOTROPIC)
log_q_ratio_barker<-function(x,y,sigma){
  beta1<-  c(-g_prime(y)*(x-y))
  beta2<-  c(-g_prime(x)*(y-x))
  return(sum(# compute acceptance with log_sum_exp trick for numerical stability
    -(pmax(beta1,0)+log1p(exp(-abs(beta1))))+
      (pmax(beta2,0)+log1p(exp(-abs(beta2))))
  ))
}
rbarker <-function(c,sigma){
  z<-sigma*rnorm(n=length(c)) # rnorm could be replaced with any symm kernel
  b<-2*(runif(n=length(c))< 1/(1+exp(-c*z)))-1
  return(z*b)
}
# the following function efficiently computes ESJD in 
# stationarity when one is able to do iid samples from the target
grad_MCMC_ESJD<-function(T=T,sigma=sigma,n=n,method=NULL){
  stopifnot(method=="MALA"|method=="RWM"|method=="Barker")
  if(method=="MALA"){
    log_q_ratio<-log_q_ratio_mala
    rprop<-rmala
  }
  if(method=="MALTA"){
    log_q_ratio<-log_q_ratio_malta
    rprop<-rmalta
  }
  if(method=="RWM"){
    log_q_ratio<-log_q_ratio_rw
    rprop<-rrw
  }
  if(method=="Barker"){
    log_q_ratio<-log_q_ratio_barker
    rprop<-rbarker
  }
  ESJD_vec<-rep(0,n)
  for (t in 1:T){
    x<-rtarget(n=n-1+1) # this function performs iid sampling from the target
    y<-x+rprop(g_prime(x),sigma)
    ap<-exp(log_f_ratio(x,y)+log_q_ratio(x,y,sigma))
    if (runif(1)<ap){
      ESJD_vec<-ESJD_vec+(y-x)^2/T      
    }
  }
  return(ESJD_vec)
}

##################################################
####### PROPOSALS WITH DIAGONAL ADAPTATION #######
##################################################
# this function implements RWM/MALA/Barker with diagonal covariance adaptation
diag_grad_MCMC<-function(T=T,n=n,method=NULL,start_x,start_sigma=NULL){
  stopifnot(method=="MALA"|method=="RWM"|method=="Barker")
  sigma<-start_sigma
  if(method=="MALA"){
    log_q_ratio_prec<-log_q_ratio_diag_mala
    rprop<-rmala_diag
    target_ap<-0.574
    if(is.null(sigma)){sigma<-2.4/sqrt(n^(1/3))}
  }
  if(method=="RWM"){
    log_q_ratio_prec<-log_q_ratio_diag_rw
    rprop<-rrw_diag
    target_ap<-0.234
    if(is.null(sigma)){sigma<-2.4/sqrt(n)}
  }
  if(method=="Barker"){
    log_q_ratio_prec<-log_q_ratio_diag_barker
    rprop<-rbarker_diag
    target_ap<-0.4
    if(is.null(sigma)){sigma<-2.4/sqrt(n^(1/3))}
  }
  # SET STARTING PARAMETERS
  x<-start_x
  t<-1
  sigma_vec<-rep(NA,T-1)
  x_samples<-matrix(NA,nrow = T,ncol = n)
  x_samples[1,]<-x
  sigma_t<-matrix(NA,nrow = T,ncol = n)
  ## params for adaption
  diag_var<-rep(1,n)
  x_means<-x
  sigma_vec[1]<-sigma
  sigma_t[1,]<-diag_var
  # RUNNING MCMC LOOP
  for (t in 2:T){
    y<-x+rprop(g_prime(x),sigma,diag_sd=sqrt(diag_var))
    if(log_f_ratio(x,y)<= -300){# if statement added for numerical stability
      ap<-0
    }else{
      ap<-min(1,exp(log_f_ratio(x,y)+log_q_ratio_prec(x,y,sigma,diag_sd=sqrt(diag_var))))
    }
    if (runif(1)<ap){
      x<-y
    }
    # store things for output
    x_samples[t,]<-x
    # do adaptation
    # 1- adapt global scale
    log_sigma_2<-log(sigma^2)+gamma(1+t)*(ap-target_ap)
    sigma<-sqrt(exp(log_sigma_2))
    # 2- adapt means
    x_means<-x_means+gamma(1+t)*(x-x_means)
    # 3- adapt diagonal covariance
    diag_var<-diag_var+gamma(1+t)*(c((x-x_means)^2)-diag_var)
    # store adaptation parameters
    sigma_vec[t]<-sigma
    sigma_t[t,]<-diag_var
  }
  return(list(x_samples=x_samples,sigma_vec=sigma_vec,sigma_t=sigma_t))
}
# function specifying the decay of the learning rate for the Robbins-Monroe algorithm in the adaptation
gamma<-function(t,kappa=0.6){
  return(t^(-kappa))
}
## MALA PROPOSAL AND ACCEPTANCE PROBABILITY (DIAGONAL)
log_q_ratio_diag_mala<-function(x,y,sigma,diag_sd){
  return(sum(
    -(x-y-(sigma*diag_sd)^2*g_prime(y)/2)^2/(2*(sigma*diag_sd)^2)+
      (y-x-(sigma*diag_sd)^2*g_prime(x)/2)^2/(2*(sigma*diag_sd)^2)
  ))}
rmala_diag <-function(c,sigma,diag_sd){
  return(
    (sigma*diag_sd)*rnorm(n=length(c))+(sigma*diag_sd)^2*c/2
  )
}
## RANDOM WALK PROPOSAL AND ACCEPTANCE PROBABILITY (DIAGONAL)
log_q_ratio_diag_rw<-function(x,y,sigma,diag_sd){return(0)}
rrw_diag<-function(c,sigma,diag_sd){return(sigma*rnorm(n=length(c),mean = 0,sd = diag_sd))}
## BARKER PROPOSAL AND ACCEPTANCE PROBABILITY (DIAGONAL)
log_q_ratio_diag_barker<-function(x,y,sigma,diag_sd){
  beta1<-  c(-g_prime(y)*(x-y))
  beta2<-  c(-g_prime(x)*(y-x))
  return(sum(# compute acceptance with log_sum_exp trick for numerical stability
    -(pmax(beta1,0)+log1p(exp(-abs(beta1))))+
      (pmax(beta2,0)+log1p(exp(-abs(beta2))))
  ))
}
rbarker_diag<-function(c,sigma,diag_sd){
  z<-sigma*rnorm(n=length(c),mean = 0,sd = diag_sd)
  b<-2*(runif(n=length(c))< 1/(1+exp(-c*z)))-1
  return(z*b)
}

###################################################################
### DEFINE THE 4 SCENARIOS USED IN HETEROGENEOUS TARGETS SECTION  #####
###################################################################
set_target<-function(scenario){
  print(paste("Set target as in scenario ",scenario))
  if(scenario==1){#Gaussian with one small scale
    n<<-100 # number of dimensions
    Sigma_targ<<-diag(rep(1,n))
    Sigma_targ[1,1]<<-0.01^2
    Prec_target<<-solve(Sigma_targ)
    log_f_ratio<<-function(x,y){return(-0.5*(
      matrix(y,nrow = 1,ncol = n)%*%Prec_target%*%matrix(y,nrow = n,ncol = 1)-
        matrix(x,nrow = 1,ncol = n)%*%Prec_target%*%matrix(x,nrow = n,ncol = 1)
    ))}
    # g_prime computes the gradient of the log_target
    g_prime<<-function(x){return(-c(matrix(x,nrow = 1,ncol = n)%*%Prec_target))}
  }
  if(scenario==2){#Gaussian with random scales
    n<<-100
    scales<<-exp(rnorm(n,mean = 0,sd = 1))
    variances<<-scales^2
    Sigma_targ<<-diag(variances)
    Prec_target<<-solve(Sigma_targ)
    log_f_ratio<<-function(x,y){return(-0.5*(
      matrix(y,nrow = 1,ncol = n)%*%Prec_target%*%matrix(y,nrow = n,ncol = 1)-
        matrix(x,nrow = 1,ncol = n)%*%Prec_target%*%matrix(x,nrow = n,ncol = 1)
    ))}
    g_prime<<-function(x){return(-c(matrix(x,nrow = 1,ncol = n)%*%Prec_target))}
  }
  if(scenario==3){#Hyperbolic with random scales
    n<<-100
    delta2<<-0.1
    scales<<-exp(rnorm(n,mean = 0,sd = 1))
    log_f_ratio<<-function(x,y){return(-sum(
      (sqrt(delta2+(y/scales)^2)-sqrt(delta2+(x/scales)^2))
    ))}
    g_prime<<-function(x){
      x_resc<- x/scales
      grad<- -x_resc/sqrt(delta2+x_resc^2)
      return(grad/scales)
    }
  }
  if(scenario==4){#Skew-normal with random scales
    n<<-100
    scales<<-exp(rnorm(n,mean = 0,sd = 1))
    alpha<<-4
    log_f_ratio<<-function(x,y){
      return(sum(
        dnorm(y/scales,0,1, log = TRUE)+pnorm(alpha*y/scales,0,1, log.p = TRUE)-
          (dnorm(x/scales,0,1, log = TRUE)+pnorm(alpha*x/scales,0,1, log.p = TRUE))
      ))
    }
    g_prime<<-function(x){
      x_resc<- x/scales
      grad<- -x_resc+
        alpha*exp(dnorm(alpha*x_resc,0,1, log = TRUE)-pnorm(alpha*x_resc,0,1, log.p = TRUE))
      return(grad/scales)
    }
  }
}

###################################################################
### DEFINE THE HIERARCHICAL POISSON REGRESSION TARGET  #####
###################################################################
set_nested_pois_target<-function(n1=10,
                                 vmu=1/(10^2), # prior precision for global means
                                 va=1/(2^2), # prior precision for random effects
                                 true_mu=5,
                                 N=NULL, # num of obs
                                 blk.ind=NULL
){
  stopifnot(length(blk.ind)==N)
  n<<-1+n1# number of parameters to sample
  par_names<<-rep(NA,n)
  par_names[1]<<-c("mu")
  ind1<<-1+1:n1
  par_names[ind1]<<-c("eta")
  ## generate data
  true_eta<<-true_mu+rnorm(n1,mean = 0,sd = 1/sqrt(va))
  blk.list<<-lapply(X = c(1:n1),FUN = function(i){which(blk.ind==i)})
  y.vec<<-c(rpois(n = N,lambda = exp(true_eta[blk.ind])))
  # target and gradient
  log_f_ratio<<-function(x,y){
    return(
      log_f(y)-log_f(x)
    )}
  log_f<<-function(x){
    mu<-x[1]
    eta<-x[ind1]
    log_lambda<-eta[blk.ind]
    return(
      -mu^2*vmu/2-sum((eta-mu)^2)*va/2+  #log prior
        sum(-exp(log_lambda)+y.vec*(log_lambda)) #log likelihoods
    )
  }
  g_prime<<-function(x,sigma_0=sigma0){
    mu<-x[1]
    eta<-x[ind1]
    lambda.vec<-exp(eta[blk.ind])
    ly_diff<-lambda.vec-y.vec
    return(c(-mu*vmu-sum(mu-eta)*va,#grad_mu
             -(eta-mu)*va-vapply(blk.list,FUN = function(ii){sum(ly_diff[ii])},FUN.VALUE = 1)#grad_eta
    ))
  }
}

##################################
#### FUNCTION TO PLOT OUTPUT #####
##################################
plot_tuning_traceplots<-function(output_MALA,output_RWM,output_Barker){
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
}
