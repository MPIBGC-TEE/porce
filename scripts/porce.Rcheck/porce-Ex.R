pkgname <- "porce"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('porce')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
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
