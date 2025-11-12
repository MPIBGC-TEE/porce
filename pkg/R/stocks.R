#' Time-dependent stocks for a linear autonomous model
#'
#' @param model an object of class lam, a linear autonomous model 
#' @param t a scalar value of the time to obtain the solution of the model
#' @param t0 a scalar value of the initial time
#' @param x0 a vector containing the size of the compartments at the initial time
#'
#' @return a vector with the stocks at the desired time
#' @export
#' @importFrom expm expm
#'
stocks <- function(t, t0, x0, model){
             if(inherits(model, 'lam') != TRUE) stop('model must be of class lam')
             matexp<-function(M,v){expm(M)%*%v}
             intmatexp<-function(x){integrate(matexp, M=model@matrix, v=model@input)}
             xt<- matexp(M=model@matrix*(t-t0), v=x0 ) #+
#                  (integrate(f=function(s){expm(model@matrix*(t-s))%*%model@input}, lower=t0, upper=t)$value)              
                    # -1 * solve(model@matrix)%*%model@input
          return(xt)
}

EM<-EmanuelModel()
stocks(model=EM, t=10, t0=0, x0=matrix(10,nrow=5))
