#' The Terrestrial Ecosystem Carbon model of Weng and Luo
#'
#' This function returns the TECO model of Weng and Luo as an object of class lam.
#' It is mostly a constructor of the model to be used for further analysis.
#'
#' @return a model of class lam
#' @export
#'
#' @examples
#' TECOModel()
TECOModel<-function(){
  U=3.370
  b=c(0.14,0.26,0.14,rep(0,5))
  mp=(1-sum(b))/3 #missing proportions of respiration. These proportions will be added back so they will be respired directly from the pools.
  bm=c(0.14+mp,0.26+mp,0.14+mp,rep(0,5))
  
  C=diag(c(0.00258,0.0000586,0.00239,0.0109,0.00095,0.0105,0.0000995,0.0000115))
  f41=0.9-(mp/2); f43=0.2-(mp/2)
  f51=0.1-(mp/2); f52=1-mp; f53=0.8-(mp/2)
  f64=0.45; f65=0.275; f67=0.42; f68=0.45
  f75=0.275; f76=0.296
  f86=0.004;f87=0.01
  
  A=matrix(c(-1, 0, 0, 0, 0, 0, 0, 0,
             0, -1, 0, 0, 0, 0, 0, 0,
             0, 0, -1, 0, 0, 0, 0, 0,
             f41, 0, f43, -1, 0, 0, 0, 0,
             f51, f52, f53, 0, -1, 0, 0, 0, 
             0, 0, 0, f64, f65, -1, f67, f68,
             0, 0, 0, 0, f75, f76, -1, 0,
             0, 0, 0, 0, 0, f86, f87, -1),
           byrow=TRUE,nrow=8,ncol=8)
  
  # Change units to MgC ha-1 yr-1
  u_TECO=U*bm*3.65
  B_TECO=(A%*%C)*365
  
  lam(input = u_TECO, matrix=B_TECO)
}


