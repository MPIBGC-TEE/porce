
#+ message=FALSE
library(porce)
library(parallel)
library(SoilR)

no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)

setwd("~/MPIBGC-TEE/porce/simulation/ATTO/")
monthly_gpp<-read.csv("inputs/ATTO_INST_monthlyGPP_2014_2019.csv")
pool_names<-c("Foliage", "Wood", "Fine roots", "Coarse roots", "Fine litter", "CWD", "Soil (0-30 cm)")
GPPmean<- mean(unlist(monthly_gpp[-1])) * 1e-06 * 1e04 *12 
GPPsd<- sd(unlist(monthly_gpp[-1])) * 1e-06 * 1e04 *12 

set.seed(1033)
atto_gpp<-rnorm(n=nrow(modpars), mean=GPPmean, sd=GPPsd/sqrt(6))

fNHZ3<-read.csv(file="inputs/NHZ3forecast.csv")
NHZ3<-data.frame(Year=c(Hua2021$NHZone3[-(937:939),1], fNHZ3$time),
                 Delta14C=c(Hua2021$NHZone3[-(937:939),2],Delta14C_from_AbsoluteFractionModern(fNHZ3$F14C)))

pars2stocks<-function(gpp, pars){
  mod<-lam(input=inputGPP(gpp), matrix=makeB(pars))
  equilibriumStock(model=mod)
}

allObjects=ls()
clusterExport(cl, varlist=c(allObjects, "lam", "inputGPP", "makeB", "equilibriumStock", "equilibriumOutflux"))

stockList<-clusterMap(cl, fun=pars2stocks, gpp=atto_gpp, pars=asplit(modpars, MARGIN = 1))
stocksATTO<-matrix(unlist(stockList), nrow=7, ncol=1000)

years<-seq(1942,2025, by=1/2) # Make the simulations shorter in time
C14_ty<-function(gpp, pars, targetYear=2024.5){
  Mod<-Model_14(t=years,inputFluxes = inputGPP(gpp),A=makeB(pars), ivList=stocksATTO[,1],
                inputFc=BoundFc(NHZ3, lag=0, format="Delta14C"), 
                initialValF=ConstFc(c(-20, -25, -20, -25, -20, -30, -50), format="Delta14C"))
  C14<-getF14(Mod)
  C14[years==targetYear,]
}


allObjects=ls()
clusterExport(cl, varlist=c(allObjects, "getF14", "ConstFc", "Model_14", "BoundFc", "inputGPP", "makeB"))
CS14<-clusterMap(cl,fun=C14_ty, gpp=atto_gpp, pars=asplit(x=modpars, MARGIN=1))
CS14_2024<-matrix(unlist(CS14), ncol=7, nrow=1000, byrow=TRUE)

mCS14<-apply(CS14_2024, MARGIN = 2, mean) # initial condition for experiment

future<-seq(2024.5, 2034.5, by=1/12)
fAtmC14<-sapply(future, function(t) -1000 )

ctrl<-Model_14(t=future,inputFluxes = inputGPP(GPPmean),A=makeB(modpars[1,]), ivList=stocksATTO[,1],
               inputFc=BoundFc(data.frame(future, fAtmC14), lag=0, format="Delta14C"), 
               initialValF=ConstFc(mCS14, format="Delta14C"))
ctrlC14<-getF14(ctrl)
ctrlC14b<-getF14C(ctrl)
ctrlC14r<-getF14R(ctrl)

library(RColorBrewer)
pal=brewer.pal(7, "Set2")

matplot(future, ctrlC14, type="l", lty=1, lwd=2, col=pal, bty="n")
legend("right", pool_names, col=pal, lty=1, bty="n")


dbl<-Model_14(t=future,inputFluxes = inputGPP(2*GPPmean),A=makeB(modpars[1,]), ivList=stocksATTO[,1],
               inputFc=BoundFc(data.frame(future, fAtmC14), lag=0, format="Delta14C"), 
               initialValF=ConstFc(mCS14, format="Delta14C"))
dblC14<-getF14(dbl)
dblC14b<-getF14C(dbl)
dblC14r<-getF14R(dbl)

plot(future, ctrlC14b, type="l", ylim=c(-1000, 100), col=pal[1], bty="n")
lines(future, ctrlC14r, col=pal[2])
lines(future, dblC14b, col=pal[3])
lines(future, dblC14r, col=pal[4])
legend("bottomleft", c("Control ecosystem radiocarbon", "Control respired radiocarbon",
                       "Double ecosystem radiocarbon", "Double respired radiocarbon"), lty=1, col=pal, bty="n")

par(mfrow=c(2,1))
matplot(future, ctrlC14, type="l", lty=1, lwd=2, col=pal, bty="n")
legend("right", pool_names, col=pal, lty=1, bty="n")
matplot(future, dblC14, type="l", lty=1, lwd=2, col=pal, bty="n")
legend("right", pool_names, col=pal, lty=1, bty="n")
par(mfrow=c(1,1))

plot(ctrlC14r, dblC14r, type="l")
abline(0,1, lty=2)
