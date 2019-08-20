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

###################################
####### DIAGONAL ADAPTATION #######
###################################
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

##########################################
####### FULL COVARIANCE ADAPTATION #######
##########################################
# this function implements RWM/MALA/Barker with full covariance adaptation
full_cov_grad_MCMC<-function(T=T,n=n,method=NULL,start_x,start_sigma=NULL,start_SIGMA=NULL,eps_reg=0){
  stopifnot(method=="MALA"|method=="MALTA"|method=="MALTAc"|method=="RWM"|method=="Barker")
  sigma<-start_sigma
  if(method=="MALA"){
    log_q_ratio_prec<-log_q_ratio_prec_mala
    rprop<-rmala_prec
    target_ap<-0.574
    if(is.null(sigma)){sigma<-2.4/sqrt(n^(1/3))}
  }
  if(method=="RWM"){
    log_q_ratio_prec<-log_q_ratio_prec_rw
    rprop<-rrw_prec
    target_ap<-0.234
    if(is.null(sigma)){sigma<-2.4/sqrt(n)}
  }
  if(method=="Barker"){
    log_q_ratio_prec<-log_q_ratio_prec_barker
    rprop<-rbarker_prec
    target_ap<-0.4
    if(is.null(sigma)){sigma<-2.4/sqrt(n^(1/3))}
  }
  # SET STARTING PARAMETERS
  x<-start_x
  t<-1
  sigma_vec<-rep(NA,T-1)
  x_samples<-matrix(NA,nrow = T,ncol = n)
  x_samples[1,]<-x
  ## parameters for preconditioning and adaption
  if(is.null(start_SIGMA)){
    SIGMA<-diag(n)
  }else{    
    SIGMA<-start_SIGMA
  }
  C<-chol(SIGMA+eps_reg*diag(n)) # if eps_reg>0 there is regularization. Default is 0
  C_inv<-solve(C)
  x_means<-x
  sigma_vec[1]<-sigma
  # RUNNING MCMC LOOP
  for (t in 2:T){
    y<-x+rprop(g_prime(x),sigma,C,C_inv)
    if(log_f_ratio(x,y)<= -300){# if statement added for numerical stability
      ap<-0
    }else{
      ap<-min(1,exp(log_f_ratio(x,y)+log_q_ratio_prec(x,y,sigma,C,C_inv)))
    }
    if (runif(1)<ap){
      x<-y
    }
    # do adaptation of global scale
    log_sigma_2<-log(sigma^2)+gamma(1+t)*(ap-target_ap)
    sigma<-sqrt(exp(log_sigma_2))
    # adapt means with robbins-monroe
    x_means<-x_means+gamma(1+t)*(x-x_means)
    # adapt covanriance with robbins-monroe
    SIGMA_0<-matrix(x-x_means,nrow=n,ncol=1)%*%matrix(x-x_means,nrow=1,ncol=n)
    SIGMA<-SIGMA+gamma(1+t)*(SIGMA_0-SIGMA)
    C<-chol(SIGMA+eps_reg*diag(n)) # if eps_reg>0 there is regularization. Default is 0
    C_inv<-solve(C)
    # store samples
    x_samples[t,]<-x
    # store adaptation parameters
    sigma_vec[t]<-sigma
  }
  return(list(x_samples=x_samples,sigma_vec=sigma_vec,SIGMA=SIGMA))
}
## MALA PROPOSAL AND ACCEPTANCE PROBABILITY (PRECONDITIONED)
log_q_ratio_prec_mala<-function(x,y,sigma,C,C_inv){
  return(
    log_q_prec_mala(y,x,sigma,C,C_inv)-log_q_prec_mala(x,y,sigma,C,C_inv)
  )
}
log_q_prec_mala<-function(x,y,sigma,C,C_inv){
  return(
    dmvnorm(x=c(y) , mean = as.vector(x+t(C)%*%C%*%g_prime(x)) , sigma = t(C)%*%C, log = TRUE
    )
  )
}
rmala_prec <-function(c,sigma,C,C_inv){
  return( t(C)%*%C%*%c + t(C) %*% rnorm(length(c)) )
}
## RANDOM WALK PROPOSAL AND ACCEPTANCE PROBABILITY (PRECONDITIONED)
log_q_ratio_prec_rw<-function(x,y,sigma,C,C_inv){return(0)}
rrw_prec<-function(c,sigma,C,C_inv){
  c_prec<-C%*%c
  return( c( t(C)%*%rrw(c_prec,sigma)  ) )
}
## BARKER PROPOSAL AND ACCEPTANCE PROBABILITY (PRECONDITIONED)
log_q_ratio_prec_barker<-function(x,y,sigma,C,C_inv){
  c_x_prec<-c(c(g_prime(x))%*%t(C))
  c_y_prec<-c(c(g_prime(y))%*%t(C))
  z<-t(C_inv)%*%(y-x)
  beta1<-  c(-c_y_prec*(-z))
  beta2<-  c(-c_x_prec*z)
  return(sum(# compute acceptance with log_sum_exp trick for numerical stability
    -(pmax(beta1,0)+log1p(exp(-abs(beta1))))+
      (pmax(beta2,0)+log1p(exp(-abs(beta2))))
  ))
}
rbarker_prec <-function(c,sigma,C,C_inv){
  c_prec<-c(c(c)%*%t(C))
  return(  t(C) %*% rbarker(c_prec,sigma)   )
}


