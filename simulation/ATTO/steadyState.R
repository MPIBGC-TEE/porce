# Steady-state carbon stocks, respiration, and radiocarbon values for ATTO
library(porce)
library(parallel)
no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)

setwd("~/MPIBGC-TEE/porce/simulation/ATTO/")

pool_names<-c("Foliage", "Wood", "Fine roots", "Coarse roots", "Fine litter", "CWD", "Soil (0-30 cm)")
monthly_gpp<-read.csv("inputs/ATTO_INST_monthlyGPP_2014_2019.csv")

matplot(as.Date(monthly_gpp[,1]), monthly_gpp[,-1]* 1e-06 * 1e04 *12, type="l", lty=1, col=1:3, ylab="GPP (MgC ha-1 yr-1)", xlab="Calendar year",
        ylim=c(0,50), bty="n")

GPPmean<- mean(monthly_gpp$GPP_U50_f) * 1e-06 * 1e04 *12 # gC/m2/month (1Mg/10e6 g) * (1e04 m2/1ha) * (12 month/1yr) -> MgC/ha/yr
GPPsd<- sd(monthly_gpp$GPP_U50_f) * 1e-06 * 1e04 *12 # gC/m2/month (1Mg/10e6 g) * (1e04 m2/1ha) * (12 month/1yr) -> MgC/ha/yr

atto_gpp<-rnorm(n=nrow(modpars), mean=GPPmean, sd=GPPsd)
hist(atto_gpp)

ATTOmodel<-lam(input=inputGPP(GPPmean), matrix=makeB(modpars[1,]))
stocksATTO<-equilibriumStock(model=ATTOmodel)

pars2stocks<-function(gpp, pars){
  mod<-lam(input=inputGPP(gpp), matrix=makeB(pars))
  equilibriumStock(model=mod)
}

pars2stocks(gpp=atto_gpp[1], pars=modpars[1,])

allObjects=ls()
clusterExport(cl, varlist=c(allObjects, "lam", "inputGPP", "makeB", "equilibriumStock"))
stockList<-clusterMap(cl, fun=pars2stocks, gpp=atto_gpp, pars=asplit(modpars, MARGIN = 1))

stocksATTO<-matrix(unlist(stockList), nrow=7, ncol=1000)

boxplot(stocksATTO, use.cols = FALSE, range=0, horizontal=TRUE, names=pool_names, 
        xlab=expression(paste("Carbon stock (Mg C h", a^-1, ")")))
