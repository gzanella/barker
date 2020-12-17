rbarker <-function(grad,sigma){
  z<-sigma*rnorm(n=length(grad)) # rnorm could be replaced with any symm kernel
  b<-2*(runif(n=length(grad))< 1/(1+exp(-grad*z)))-1
  return(z*b)
}
