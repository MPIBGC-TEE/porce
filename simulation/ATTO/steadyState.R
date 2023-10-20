# Steady-state carbon stocks, respiration, and radiocarbon values for ATTO
library(porce)
library(parallel)
no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)

setwd("~/MPIBGC-TEE/porce/simulation/ATTO/")

pool_names<-c("Foliage", "Wood", "Fine roots", "Coarse roots", "Fine litter", "CWD", "Soil (0-30 cm)")
monthly_gpp<-read.csv("inputs/ATTO_INST_monthlyGPP_2014_2019.csv")

matplot(as.Date(monthly_gpp[,1]), monthly_gpp[,-1]* 1e-06 * 1e04 *12, type="l", 
        lty=1, col=1:3, ylab=expression(paste("GPP (MgC h", a^-1, " y", r^-1,")")), xlab="Calendar year",
        ylim=c(0,50), bty="n")

GPPmean<- mean(unlist(monthly_gpp[-1])) * 1e-06 * 1e04 *12 # gC/m2/month (1Mg/10e6 g) * (1e04 m2/1ha) * (12 month/1yr) -> MgC/ha/yr
GPPsd<- sd(unlist(monthly_gpp[-1])) * 1e-06 * 1e04 *12 # gC/m2/month (1Mg/10e6 g) * (1e04 m2/1ha) * (12 month/1yr) -> MgC/ha/yr

atto_gpp<-rnorm(n=nrow(modpars), mean=GPPmean, sd=GPPsd/6)
hist(atto_gpp)

ATTOmodel<-lam(input=inputGPP(GPPmean), matrix=makeB(modpars[1,]))
stocksATTO<-equilibriumStock(model=ATTOmodel)

pars2stocks<-function(gpp, pars){
  mod<-lam(input=inputGPP(gpp), matrix=makeB(pars))
  equilibriumStock(model=mod)
}

pars2stocks(gpp=atto_gpp[1], pars=modpars[1,])

allObjects=ls()
clusterExport(cl, varlist=c(allObjects, "lam", "inputGPP", "makeB", "equilibriumStock", "equilibriumOutflux"))
stockList<-clusterMap(cl, fun=pars2stocks, gpp=atto_gpp, pars=asplit(modpars, MARGIN = 1))

stocksATTO<-matrix(unlist(stockList), nrow=7, ncol=1000)

boxplot(stocksATTO, use.cols = FALSE, range=0, horizontal=TRUE, names=pool_names, 
        xlab=expression(paste("Carbon stock (Mg C h", a^-1, ")")))

pars2resp<-function(gpp, pars){
  mod<-lam(input=inputGPP(gpp), matrix=makeB(pars))
  equilibriumOutflux(model=mod)
}

respList<-clusterMap(cl, fun=pars2resp, gpp=atto_gpp, pars=asplit(modpars, MARGIN = 1))
respATTO<-matrix(unlist(respList), nrow=7, ncol=1000)

boxplot(respATTO, use.cols = FALSE, range=0, horizontal=TRUE, names=pool_names, 
        xlab=expression(paste("Respiration (Mg C h", a^-1, " y", r^-1, ")")))

abvg_pools<-match(c("Foliage", "Wood"), table=pool_names)
blwg_pools<-match(c("Fine roots", "Coarse roots"), table=pool_names)
necromass_pools<-match(c("Fine litter", "CWD"), table=pool_names)

abvg<-colSums(stocksATTO[abvg_pools,])
blwg<-colSums(stocksATTO[blwg_pools,])
necro<-colSums(stocksATTO[necromass_pools,])

boxplot(list(abvg, blwg, necro), names=c("Aboveground C", "Belowground C", "Necromass C"), ylim=c(0,300), range=0,
        ylab=expression(paste("Respiration (Mg C h", a^-1, ")")))

toc<-colSums(stocksATTO)

autotroph_pools=match(c("Foliage", "Wood", "Fine roots", "Coarse roots"), table=pool_names)
heterotroph_pools=match(c("Fine litter", "CWD", "Soil (0-30 cm)"), table=pool_names)

Ra<-colSums(respATTO[autotroph_pools,])
Rh<-colSums(respATTO[heterotroph_pools,])
Re<-colSums(respATTO)

boxplot(list(Ra,Rh,Re), names=c("Autotrophic", "Heterotrophic", "Ecosystem"), ylim=c(0,60), range=0,
        ylab=expression(paste("Respiration (Mg C h", a^-1, " y", r^-1, ")")))


npp<-atto_gpp-Ra
NG<-npp/atto_gpp

par(mfrow=c(2,1), mar=c(4,4,1,1))
hist(npp, main="", xlab=expression(paste("Net primiary production NPP (Mg C h", a^-1, " y", r^-1, ")")))
hist(NG, main="", xlab="NPP:GPP")
par(mfrow=c(1,1))


image(x=1:50,y=1:20,z=matrix(toc,50,20), xlab="x coordinate", ylab="y coordinate")

library(SoilR)

fNHZ3<-read.csv(file="inputs/NHZ3forecast.csv")
NHZ3<-data.frame(Year=c(Hua2021$NHZone3[-(937:939),1], fNHZ3$time),
                 Delta14C=c(Hua2021$NHZone3[-(937:939),2],Delta14C_from_AbsoluteFractionModern(fNHZ3$F14C)))

years<-seq(1942,2029, by=1/12)
C14Model<-Model_14(t=years,inputFluxes = inputGPP(GPPmean),A=makeB(modpars[1,]), ivList=stocksATTO[,1],
                   inputFc=BoundFc(NHZ3, lag=0, format="Delta14C"), 
                   initialValF=ConstFc(c(-20, -25, -20, -25, -20, -30, -50), format="Delta14C"))
C14t<-getF14(C14Model)

plot(NHZ3, type="l", col=1, lwd=2)
matlines(years, C14t, lty=1, col=rainbow(n=7))
legend("topright", c("Atmosphere",pool_names), col=c(1, rainbow(n=7)), lty=1, lwd=c(2,rep(1,7)), bty="n")


R14t<-getF14R(C14Model)

plot(NHZ3, type="l", col=1, lwd=2, ylab=expression(paste(Delta^14, "C (per mil)")), bty="n")
lines(years,R14t, col=2)
legend("topright", c("Atmosphere", "Respired CO2"), lty=1, col=1:2, bty="n")


R14_ty<-function(gpp, pars, targetYear=2019.5){
  Mod<-Model_14(t=years,inputFluxes = inputGPP(gpp),A=makeB(pars), ivList=stocksATTO[,1],
                inputFc=BoundFc(NHZ3, lag=0, format="Delta14C"), 
                initialValF=ConstFc(c(-20, -25, -20, -25, -20, -30, -50), format="Delta14C"))
  R14<-getF14R(Mod)
  R14[years==targetYear]
}


R14_ty(gpp=atto_gpp[1], pars=modpars[1,])

allObjects=ls()
clusterExport(cl, varlist=c(allObjects, "getF14R", "ConstFc", "Model_14", "BoundFc", "inputGPP", "makeB"))

ER14<-mapply(FUN=R14_ty, gpp=atto_gpp[1:3], pars=asplit(x=modpars[1:3,], MARGIN=1))

ER14<-clusterMap(cl,fun=R14_ty, gpp=atto_gpp, pars=asplit(x=modpars, MARGIN=1))

ER14_2019<-unlist(ER14)

pdf("Expected14CO2.pdf", encoding = 'WinAnsi.enc')
hist(ER14_2019, freq=FALSE, main="", xlab=expression(paste(Delta^14, "C (\u2030)")))
dev.off()

mean(ER14_2019)
median(ER14_2019)
range(ER14_2019)


C14_ty<-function(gpp, pars, targetYear=2019.5){
  Mod<-Model_14(t=years,inputFluxes = inputGPP(gpp),A=makeB(pars), ivList=stocksATTO[,1],
                inputFc=BoundFc(NHZ3, lag=0, format="Delta14C"), 
                initialValF=ConstFc(c(-20, -25, -20, -25, -20, -30, -50), format="Delta14C"))
  C14<-getF14(Mod)
  C14[years==targetYear]
}

clusterExport(cl, varlist=c(allObjects, "getF14R", "ConstFc", "Model_14", "BoundFc", "inputGPP", "makeB", "getF14"))
RS<-clusterMap(cl,fun=C14_ty, gpp=atto_gpp, pars=asplit(x=modpars, MARGIN=1))
radiocarbonStocks<-matrix(unlist(RS), ncol=7, nrow=1000, byrow = TRUE)
boxplot(radiocarbonStocks, names=pool_names, range=0, ylab=expression(paste(Delta^14, "C (\u2030)")))
