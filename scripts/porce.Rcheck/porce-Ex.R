pkgname <- "porce"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('porce')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("EmanuelModel")
### * EmanuelModel

flush(stderr()); flush(stdout())

### Name: EmanuelModel
### Title: Terrestrial carbon model of Emanuel
### Aliases: EmanuelModel

### ** Examples

EmanuelModel()



cleanEx()
nameEx("convolutionfun")
### * convolutionfun

flush(stderr()); flush(stdout())

### Name: convolutionfun
### Title: Convolution function
### Aliases: convolutionfun

### ** Examples

tms<-seq(0,10, by=0.1)
cosconvfun<-convolutionfun(f=cos, g=cos) # convolve cosine function with itself
solfun<-function(t){ (t*cos(t) + sin(t))/2} # analytical solution (Braun 1993, Dif Eq and App, Springer, p. 254)
x1<-sapply(X=tms, FUN=cosconvfun)
x2<-sapply(X=tms, FUN=solfun)

plot(tms, x1, type="l")
lines(tms, x2, col=2)



cleanEx()
nameEx("erf")
### * erf

flush(stderr()); flush(stdout())

### Name: erf
### Title: Equilibrium response function
### Aliases: erf

### ** Examples

EmanuelERF<-erf(model=EmanuelModel())
tms<-seq(0,100)
em<-sapply(tms, EmanuelERF)
plot(tms, em, type="l")




cleanEx()
nameEx("inputGPP")
### * inputGPP

flush(stderr()); flush(stdout())

### Name: inputGPP
### Title: Input vector from a scalar gpp value
### Aliases: inputGPP

### ** Examples

inputGPP(25, 7)



cleanEx()
nameEx("irf")
### * irf

flush(stderr()); flush(stdout())

### Name: irf
### Title: Impulse response function
### Aliases: irf

### ** Examples

EmanuelIRF<-irf(model=EmanuelModel())
tms<-seq(0,100)
em<-sapply(tms, EmanuelIRF)
plot(tms, em, type="l")




cleanEx()
nameEx("lam-class")
### * lam-class

flush(stderr()); flush(stdout())

### Name: lam-class
### Title: Linear autonomous model
### Aliases: lam-class lam

### ** Examples

toyModel<-lam(input=c(1,2,3), matrix=diag(-1,3,3))



cleanEx()
nameEx("makeB")
### * makeB

flush(stderr()); flush(stdout())

### Name: makeB
### Title: Compartmental matrix from a set of prior parameters of Porce
###   model
### Aliases: makeB

### ** Examples

makeB(pars=modpars[1,])



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
