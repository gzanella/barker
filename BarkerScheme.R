## this is a simple function to illustrate how to sample from the Barker proposal (without preconditioning)
## grad is the target log-posterior gradient (vector)
## sigma is the stepsize (scalar)
rbarker <-function(grad,sigma){
  z<-rnorm(n=length(grad),mean = 0,sd = sigma) 
  b<-2*(runif(n=length(grad))< 1/(1+exp(-grad*z)))-1
  return(z*b)
}

## A SPEEDED-UP VERSION
## above we give a N(0,sigma^2) distribution to the noize z. This can be replaced with basically any other kernel.
## A speed-up of the sampler's mixing (roughly by a factor of 2.4 in high dimensions) is obtained by using 
## a noise z with distribution centered around +/-sigma rather than at 0. In the following code we use a N(sigma,(0.1*sigma)^2)
## noise rather than a N(0,sigma^2) one. We suggest using this (or something similar) in practice.
rbarker <-function(grad,sigma){
  z<-rnorm(n=length(grad),mean = sigma,sd = 0.1*sigma)
  b<-2*(runif(n=length(grad))< 1/(1+exp(-grad*z)))-1
  return(z*b)
}
