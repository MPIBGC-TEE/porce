library(porce)
library(forecast)
library(SoilR)

setwd("~/MPIBGC-TEE/porce/simulation/ATTO/")
monthly_gpp<-read.csv("inputs/ATTO_INST_monthlyGPP_2014_2019.csv")
gpp_ts<-lapply(X=as.list(monthly_gpp[-1]), FUN=ts, start=c(2014,1), frequency = 12)
gpp_ets<-lapply(X=gpp_ts, FUN=ets)
pool_names<-c("Foliage", "Wood", "Fine roots", "Coarse roots", "Fine litter", "CWD", "Soil (0-30 cm)")

lapply(gpp_ets, plot)
ny<-10; freq<-12

gpp_forecast<-lapply(gpp_ets, FUN=forecast, h=ny*freq)
gpp_sim<-sapply(gpp_forecast, function(x) x$mean)

years<-seq(1/freq, ny, by=1/freq)
matplot(years, gpp_sim, type="l", lty=1, col=1:3, bty="n")

xss=equilibriumStock(lam(input=inputGPP(mean(gpp_sim[,1])), matrix=makeB(modpars[1,])))

influxFunc<-function(t, series=data.frame(years,gpp_sim[,1])){
  infunc<-splinefun(x=series[,1],y=series[,2])
  inputGPP(gpp=infunc(t))
}

mod<-Model(t=years,A=makeB(modpars[1,]),ivList = as.numeric(xss), inputFluxes = influxFunc)
Ct<-getC(mod)
Rt<-getReleaseFlux(mod)

autotroph_pools=match(c("Foliage", "Wood", "Fine roots", "Coarse roots"), table=pool_names)
heterotroph_pools=match(c("Fine litter", "CWD", "Soil (0-30 cm)"), table=pool_names)

plot(years, rowSums(Rt), type="l", ylab=expression(paste("Flux (Mg C h", a^-1, " y", r^-1, ")")), ylim=c(0,400), col=2, bty="n")
lines(years, gpp_sim[,1], col=3)
lines(years, rowSums(Rt[,autotroph_pools]),col=4)
lines(years, rowSums(Rt[,heterotroph_pools]),col=5)
legend("topright", c("Gross primary production", "Ecosystem respiration",
                     "Autotrophic respiration", "Heterotrophic respiration"), lty=1, col=c(3,2,4,5),bty="n")

library(parallel)
no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)

pars2resp<-function(pars){
  mod<-Model(t=years,A=makeB(pars),ivList = as.numeric(xss), inputFluxes = influxFunc)
  getReleaseFlux(mod)
}
allObjects=ls()
clusterExport(cl, varlist=c(allObjects, "lam", "inputGPP", "makeB", "Model", "getReleaseFlux"))
respList<-clusterMap(cl, fun=pars2resp, pars=asplit(modpars, MARGIN = 1))

EcoResp<-sapply(respList, FUN=function(x) rowSums(x))
meanER<-apply(X=EcoResp, MARGIN=1, FUN=mean)
sdER<-apply(X=EcoResp, MARGIN=1, FUN=sd)

AutResp<-sapply(respList, FUN=function(x) rowSums(x[,autotroph_pools]))
meanAR<-apply(X=AutResp, MARGIN=1, FUN=mean)
sdAR<-apply(X=AutResp, MARGIN=1, FUN=sd)

HetResp<-sapply(respList, FUN=function(x) rowSums(x[,heterotroph_pools]))
meanHR<-apply(X=HetResp, MARGIN=1, FUN=mean)
sdHR<-apply(X=HetResp, MARGIN=1, FUN=sd)

plot(years, gpp_sim[,1], type="l", ylab=expression(paste("Flux (Mg C h", a^-1, " y", r^-1, ")")), ylim=c(0,400), col=3, bty="n")
polygon(x=c(years,rev(years)), y=c(meanER+sdER, rev(meanER-sdER)), col=2, border=NA)
polygon(x=c(years,rev(years)), y=c(meanAR+sdAR, rev(meanAR-sdAR)), col=4, border=NA)
polygon(x=c(years,rev(years)), y=c(meanHR+sdHR, rev(meanHR-sdHR)), col=5, border=NA)
lines(years, gpp_sim[,1], col=3)
legend("topright", c("Gross primary production", "Ecosystem respiration",
                     "Autotrophic respiration", "Heterotrophic respiration"), lty=1, col=c(3,2,4,5),bty="n")

pal<-rainbow(7)
matplot(years,Rt, type="l", lty=1, ylim=c(0,200), col=pal, bty="n")
legend("topright", legend = pool_names, col=pal, lty=1, bty="n")
