#' Equilibrium stocks for a linear autonomous model
#'
#' @param model an object of class lam, a linear autonomous model 
#'
#' @return a vector with the equilibrium stocks for all compartments
#' @export
#'
equilibriumStock <- function(model){
                      if(inherits(model, 'lam') != TRUE) stop('model must be of class lam')
                     -1 * solve(model@matrix)%*%model@input
}

