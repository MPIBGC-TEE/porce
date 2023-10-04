#' makeB makes a compartmental matrix B from a vector of model parameters
#' 
#' @param pars A vector of named model parameters
#' 
makeB=function(pars){
  B=diag(-1*pars[1:7])
  B[2,1]=pars["alpha21"]*pars["k1"]
  B[3,1]=pars["alpha31"]*pars["k1"]
  B[4,1]=pars["alpha41"]*pars["k1"]
  B[5,1]=pars["alpha51"]*pars["k1"]
  B[5,3]=pars["alpha53"]*pars["k3"]
  B[6,2]=pars["alpha62"]*pars["k2"]
  B[6,4]=pars["alpha64"]*pars["k4"]
  B[7,5]=pars["alpha75"]*pars["k5"]
  B[7,6]=pars["alpha76"]*pars["k6"]
  return(B)
}
