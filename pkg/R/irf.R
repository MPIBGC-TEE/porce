#' Impulse response function
#'
#' Builds an impulse response function from a model object of class lam
#'
#' @param model A model of class lam
#'
#' @return A function that takes time as main argument
#' @export
#' @importFrom expm expm
#'
#' @examples
#' EmanuelIRF<-irf(model=EmanuelModel())
#' tms<-seq(0,100)
#' em<-sapply(tms, EmanuelIRF)
#' plot(tms, em, type="l")
#' 
irf<-function(model){
  u<-model@input/sum(model@input)
  B<-model@matrix
  z<- -1 * colSums(B)
  function(t) as.numeric(z%*%expm((t*B))%*%u)
}
