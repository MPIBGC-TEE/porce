#' Steady state C release (respiration) from the system
#' 
#' This function solves for the steady-state C release for each pool by given values of GPP 
#' and parameters of a compartmental matrix
#' 
#' @param gpp Gross primary production. A scalar value
#' @param pars A vector model parameters from which a compartmental matrix can be built
#' 
steady_release<-function(gpp, pars){
  B<-makeB(pars)
  u<-matrix(c(gpp, rep(0,(ncol(B)-1))),ncol=1)
  x<- -1*solve(B)%*%u
  r<- -1*colSums(B)
  r*x
}
