#' Compartmental matrix from a vector of named parameters
#' 
#' This function builds a compartmental matrix using a vector of
#' parameter values. The parameter must be named using a name of the form "ki" 
#' for rates for compartment i, and transfer coefficients "alphaij" with j the
#' number of the donor compartment and i the number of the receiver compartment.
#'
#' @param pars a numeric vector of appropriately named parameter values
#'
#' @return A compartmental matrix
#' @export
#'
#' @examples
#' names(modpars[1,])
#' makeB(pars=modpars[1,])
makeB=function(pars){
  # Old version for a 7 pool model with specific structure
  # B=diag(-1*pars[1:7])
  # B[2,1]=pars["alpha21"]*pars["k1"]
  # B[3,1]=pars["alpha31"]*pars["k1"]
  # B[4,1]=pars["alpha41"]*pars["k1"]
  # B[5,1]=pars["alpha51"]*pars["k1"]
  # B[5,3]=pars["alpha53"]*pars["k3"]
  # B[6,2]=pars["alpha62"]*pars["k2"]
  # B[6,4]=pars["alpha64"]*pars["k4"]
  # B[7,5]=pars["alpha75"]*pars["k5"]
  # B[7,6]=pars["alpha76"]*pars["k6"]
  # return(B)
  
  # New implementation for any dimension with any structure
  i_k<-grep("k",x=names(pars))
  m<-length(i_k)
  k_pos<-as.numeric(substr(names(pars[i_k]), 2, length(pars)))
  ks<-pars[order(k_pos)]
  
  K<-diag(abs(ks), ncol=m, nrow = m)
  
  i_alpha<-grep("alpha",x=names(pars))
  row_alpha<-as.numeric(substr(names(pars[i_alpha]), 6, 6))
  col_alpha<-as.numeric(substr(names(pars[i_alpha]), 7, 7))
  
  T<-diag(-1,nrow=m, ncol=m)
  
  for(i in 1:length(i_alpha)){
    T[row_alpha[i], col_alpha[i]]<-pars[i_alpha[i]]
  }
  return(T%*%K)
  
}
