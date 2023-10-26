#' Terrestrial carbon model of Emanuel
#'
#' This function returns the model of Emanuel as an object of class lam.
#' It is mostly a constructor of the model to be used for further analysis.
#'
#' @return a model of class lam
#' @export
#'
#' @examples
#' EmanuelModel()
EmanuelModel<-function(){
  x=matrix(c(37, 452, 69, 81, 1121),5,5, byrow=TRUE)
  Fl=matrix(c(-(21+31+25), 0, 0, 0, 0,
            31, -(14+15+2), 0, 0, 0,
            0, 0, -(18+12+6),0, 0,
            21, 15, 12, -(45+3), 0,
            0, 2, 6, 3, -(11+0.5)), 5,5, byrow=TRUE)
  B_emanuel<-Fl/x
  u_emanuel<-c(77, 0, 36, 0, 0)

  lam(input = u_emanuel, matrix=B_emanuel)
}


