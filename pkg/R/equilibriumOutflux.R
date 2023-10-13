#' Equilibrium output flux for a linear autonomous model
#'
#' @param model an object of class lam, a linear autonomous model 
#'
#' @return a vector with the output fluxes for all compartments
#' @export
#'
equilibriumOutflux <- function(model){
                      if(inherits(model, 'lam') != TRUE) stop('model must be of class lam')
                       x<- -1 * solve(model@matrix)%*%model@input
                       z<- -1 * colSums(model@matrix)
                       z*x
}

