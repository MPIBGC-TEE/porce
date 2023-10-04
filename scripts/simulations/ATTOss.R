# Simulating C stocks and radiocarbon for ATTO at steady-state

library(SoilR)
library(parallel)
no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)

pool_names<-c("Foliage", "Wood", "Fine roots", "Coarse roots", "Fine litter", "CWD", "Soil (0-30 cm)")

sapply(X=c("~/MPIBGC-TEE/porce/scripts/initial_functions/makeB.R",
          "~/MPIBGC-TEE/porce/scripts/initial_functions/steady_stocks.R",
          "~/MPIBGC-TEE/porce/scripts/initial_functions/steady_release.R"),
      FUN=source)
load("~/SOIL-R/Manuscripts/RadiocarbonPorce/Data/modpars.Rdata") # Model parameters from Sierra et al. (2021, J Ecology)

monthly_gpp<-read.csv("~/MPIBGC-TEE/porce/scripts/simulations/ATTO_INST_monthlyGPP_2014_2019.csv")

matplot(as.Date(monthly_gpp[,1]), monthly_gpp[,-1]* 1e-06 * 1e04 *12, type="l", lty=1, col=1:3, ylab="GPP (MgC ha-1 yr-1)", xlab="Calendar year",
        ylim=c(0,50), bty="n")

GPPmean<- mean(monthly_gpp$GPP_U50_f) * 1e-06 * 1e04 *12 # gC/m2/month (1Mg/10e6 g) * (1e04 m2/1ha) * (12 month/1yr) -> MgC/ha/yr
GPPsd<- sd(monthly_gpp$GPP_U50_f) * 1e-06 * 1e04 *12 # gC/m2/month (1Mg/10e6 g) * (1e04 m2/1ha) * (12 month/1yr) -> MgC/ha/yr

atto_gpp<-rnorm(n=nrow(modpars), mean=GPPmean, sd=GPPsd)
hist(atto_gpp)

#stocksATTO<-steady_stocks(gpp=GPPmean, pars=modpars[1,])

allObjects=ls()
clusterExport(cl, varlist=c(allObjects))

stocksATTO<-clusterMap(cl, fun=steady_stocks, gpp=atto_gpp, pars=asplit(modpars, MARGIN = 1))

up<-matrix(unlist(stocksATTO), nrow=7, ncol=1000)

setwd("~/MPIBGC-TEE/porce/scripts/simulations/")

pdf("stocksATTOss.pdf", height=7*sqrt(2))
boxplot(up, use.cols = FALSE, range=0, horizontal=TRUE, names=pool_names, 
        xlab=expression(paste("Carbon stock (Mg C h", a^-1, ")")))
dev.off()

respirationATTO<-clusterMap(cl, fun=steady_release, gpp=atto_gpp, pars=asplit(modpars, MARGIN = 1))
unlist_resp<-matrix(unlist(respirationATTO), nrow=7, ncol=1000)

pdf("poolRespATTOss.pdf", height=7*sqrt(2))
boxplot(unlist_resp, use.cols = FALSE, range=0, horizontal=TRUE, names=pool_names, 
        xlab=expression(paste("Respiration (Mg C h", a^-1, " y", r^-1, ")")))
dev.off()

autotroph_pools=1:4
heterotroph_pools=5:7

Ra<-colSums(unlist_resp[autotroph_pools,])
Rh<-colSums(unlist_resp[heterotroph_pools,])
Re<-colSums(unlist_resp)

pdf("respATTOss.pdf")
par(mar=c(2,4.5,1,1))
boxplot(list(Ra,Rh,Re), names=c("Ra", "Rh", "Re"), ylim=c(0,60), range=0,
        ylab=expression(paste("Respiration (Mg C h", a^-1, " y", r^-1, ")")))
dev.off()

nee<-atto_gpp-Re
boxplot(nee)
abline(h=0)

## Radiocarbon
fNHZ3<-read.csv(file="~/SOIL-R/Data/forecast14CO2/NHZ3forecast.csv")
NHZ3<-data.frame(Year=c(Hua2021$NHZone3[-(937:939),1], fNHZ3$time),
                 Delta14C=c(Hua2021$NHZone3[-(937:939),2],Delta14C_from_AbsoluteFractionModern(fNHZ3$F14C)))

## Plot
library(imager)
atto_gpp<-rnorm(n=10000, mean=GPPmean, sd=GPPsd)

graygpp<-as.cimg(matrix(atto_gpp/max(atto_gpp),100,100))
plot(graygpp)
save.image(graygpp, "atto_gpp.jpg")
