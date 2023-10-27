#' ---
#' title: "Computing ecosystem respiration using response functions"
#' author: Carlos A. Sierra
#' output: html_document
#' ---
#' 
#' Response functions can be used to compute the output signal produced as
#' a result of an input signal passing through a system. `porce` includes
#' functions to compute the impulse response function `irf` and the equilibrium
#' response function `erf` of a model of class `lam`. Here we will use this 
#' functions to show how to compute respiration fluxes for the model fo Emanuel.
#' 
#' First we prepare the packages needed for this example and a color palette. 
#' 
#+ message=FALSE
library(porce)
library(SoilR)
library(expm)
library(RColorBrewer)
pal=brewer.pal(5, "Set2")

#' The Emanuel model is already available in `porce`, we only need to store it
#' in an object and compute it's equilibrium stock.
EM<-EmanuelModel()
xss<-equilibriumStock(EM)

#' We will now create a function to perturb the values of GPP and increase it 
#' exponentially at a rate of 0.01 per year. The function `ut` returns the vector
#' of inputs to the system for every simulation year, and the object `gpp` contains
#' the numberical values
yr<-seq(0,100,by=1/12)
ut<-function(t){exp(0.01*t) * EM@input}
gpp<-sapply(yr, FUN=function(t) sum(ut(t)))

#' The `irf` and `erf` functions return the impulse response and equilibrium response
#' functions using as mail argument a model of class `lam`. The numerical evaluations
#' of the functions can be obtained by the vectorized approach of R (apply family of functions),
#' which is a convenient approach for fast computations. 
irfEM<-sapply(yr, irf(EM))
erfEM<-sapply(yr, erf(EM))

#' Although `irf` and `erf` produced the desired output, every call to these functions
#' implies a number of computations on matrices (matrix exponential and multiplication),
#' which are computationally expensive. A much simplier approach for the convolution
#' of functions is to use instead a spline function of the numerical output. 
#' In the code below we use the spline functions of the gpp and irf to obtaine the 
#' respiration response of new inputs. 

recRe<-convolutionfun(t=yr,f=splinefun(yr, gpp), g=splinefun(yr, irfEM))
ReNew<-sapply(yr,recRe)

#' The existing stocks at $t=0$ are also respired, and we use the erf 
#' to obtain the release of C from these existing stocks.

ReOld<-sum(xss)*erfEM

#' To evaluate the accuracy of these computations, we can also run a `SoilR` model
#' with time dependent inputs and compute the respiration flux. In `SoilR` the models
#' are solved using numerical solvers for a system of ordinary differential equations (ODE).
#' This approach is also a little more computational expensive, but provides a highly
#' accurate result. 

EM_SoilR<-Model(t=yr, A=EM@matrix, ivList = c(xss), inputFluxes = ut)
Rt<-getReleaseFlux(EM_SoilR)

#' Now we can visually compare the results from both approaches. We can see
#' that the response function approach yields results that are indistinguishable 
#' from the traditional approach of solving systems of ODEs. 

#+ fig.width=10, fig.height=8
par(mar=c(4,4.5,1,1), lwd=2)
plot(yr, gpp, type="l", ylim=c(0,300), xlim=c(0,100), xlab="Year",
     ylab=expression(paste("Carbon flux (Pg C ",yr^-1,")")), col=pal[1],bty="n")
lines(yr,rowSums(Rt),col=pal[2])
lines(yr,ReNew,col=pal[3])
lines(yr, ReOld, col=pal[4])
lines(yr, ReNew+ReOld, col=1, lty=2)
legend("topleft", c("GPP", "Re", "Impulse response convolution", "Equilibrium response function", "Reconstruction"),
       lty=c(rep(1,4),2),col=c(pal[-5],1),bty="n")

#' ## Transit time and age distributions
#' The impulse response function and the transit time distribution of a linear
#' autonomous model are equivalent. Similarly, the equilibrium response function
#' and the system age distribution are equivalent. Therefore, the same computations
#' above can be replaced by computations based on transit time and age distributions.
#' These distributions can be easily obtained using `SoilR` functions `transitTime`
#' and `systemAge`. 
#' 
#' One important implication of this equivalence is that the respiration response
#' of an ecosystem to a change in inputs can be obtained from the age of carbon 
#' in the respiration flux (transit time) and the age of carbon stored in an 
#' ecosystem. 
#' 