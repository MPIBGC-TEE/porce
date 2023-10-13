#' Linear autonomous model
#'
#' @slot input numeric vector with inputs for each compartment.
#' @slot matrix a compartmental matrix with dimension equal to length of input.
#'
#' @return An object of lam class
#' @export
#'
#' @examples
#' toyModel<-new("lam", input=c(1,2,3), matrix=diag(-1,3,3))

lam<-setClass(Class="lam", slots=c(input="numeric", matrix="array"))
