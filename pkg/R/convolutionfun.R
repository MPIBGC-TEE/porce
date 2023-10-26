#' Convolution function
#' 
#' Computes the convolution function between two functions f and g. 
#' Instead of returning a numerical value, this function returns a function.
#'
#' @param t a vector of time points
#' @param t0 initial time
#' @param f a function to convolve
#' @param g a function to convolve
#'
#' @return a convolution function between f and g
#' @export
#' @importFrom stats integrate
#'
#' @examples
#' tms<-seq(0,10, by=0.1)
#' cosconvfun<-convolutionfun(f=cos, g=cos) # convolve cosine function with itself
#' solfun<-function(t){ (t*cos(t) + sin(t))/2} # analytical solution (Braun 1993, Dif Eq and App, Springer, p. 254)
#' x1<-sapply(X=tms, FUN=cosconvfun)
#' x2<-sapply(X=tms, FUN=solfun)
#'
#' plot(tms, x1, type="l")
#' lines(tms, x2, col=2)

convolutionfun=function(t, t0=0, f, g){
  function(t){integrate(function(u,t){f(t-u)*g(u)}, lower=t0, upper=t, t)$value}
}

