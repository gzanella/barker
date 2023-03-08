## this is a simple function to illustrate how to sample from the Barker proposal (isotropic version, i.e. no preconditioning)
rbarker <-function(x,grad,sigma){
  # x: current location (vector)
  # grad: target log-posterior gradient (vector)
  # sigma: proposal stepsize (scalar)
  z<-rnorm(n=length(grad),mean = 0,sd = sigma) 
  b<-2*(runif(n=length(grad))< 1/(1+exp(-grad*z)))-1
  return(x+z*b)
}

## the contribution to the log-acceptance rate coming from
## the proposal, i.e. log(q(y,x)/q(x,y)), is computed as follows
log_q_ratio_barker<-function(x,y,grad_x,grad_y){
  # x: current location (vector)
  # y: proposed location (vector)
  # grad_x: target log-posterior gradient at x (vector)
  # grad_y: target log-posterior gradient at y (vector)
  beta1<-  c(-grad_y*(x-y))
  beta2<-  c(-grad_x*(y-x))
  return(sum(# compute acceptance with log_sum_exp trick for numerical stability
    -(pmax(beta1,0)+log1p(exp(-abs(beta1))))+
      (pmax(beta2,0)+log1p(exp(-abs(beta2))))
  ))
}

## A SIMPLE SPEED-UP TRICK
## Above we give a N(0,sigma^2) distribution to the noize z. This can be replaced with any other kernel (without changing the value of the acceptance function!).
## A speed-up of the sampler's mixing (roughly by a factor of 2.4 in high dimensions) is obtained by using 
## a noise z with distribution centered around +/-sigma rather than at 0. In the following code we use a N(sigma,(0.1*sigma)^2)
## noise rather than a N(0,sigma^2) one. We suggest using this (or something similar) in practice.
rbarker <-function(x,grad,sigma){
  # x: current location (vector)
  # grad: target log-posterior gradient (vector)
  # sigma: proposal stepsize (scalar)
  z<-rnorm(n=length(grad),mean = sigma,sd = 0.1*sigma)
  b<-2*(runif(n=length(grad))< 1/(1+exp(-grad*z)))-1
  return(x+z*b)
}
