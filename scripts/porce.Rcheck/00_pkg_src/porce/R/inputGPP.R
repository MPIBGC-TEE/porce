#' Input vector from a scalar gpp value
#'
#' @param gpp a scalar value of gross primary production
#' @param npools integer. Number of pools in the system
#'
#' @return a vector of npool elements with GPP as first argument
#' @export
#'
#' @examples
#' inputGPP(25, 7)
inputGPP <- function(gpp, npools=7){
  c(gpp, rep(0, npools-1))
}