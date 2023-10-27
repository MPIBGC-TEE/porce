#' Equilibrium response function
#'
#' Builds an equilibrium response function from a model object of class lam
#'
#' @param model A model of class lam
#'
#' @return A function that takes time as main argument
#' @export
#' @importFrom expm expm
#'
#' @examples
#' EmanuelERF<-erf(model=EmanuelModel())
#' tms<-seq(0,100)
#' em<-sapply(tms, EmanuelERF)
#' plot(tms, em, type="l")
#' 
erf<-function(model){
  u<-model@input
  B<-model@matrix
  z<- -1 * colSums(B)
  x<- -1*solve(B)%*%u
  function(t) as.numeric(z%*%expm((t*B))%*%(x/sum(x)))
}
