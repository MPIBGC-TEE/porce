#' ---
#' title: "Radiocarbon response to a CO~2~ fertilization experiment with fossil carbon"
#' author: Carlos A. Sierra
#' date: "2023-10-31"
#' output: html_document
#' ---
#'
#' This document contains a set of simulations to predict the radiocarbon response
#' of a tropical ecosystem, representative of the central Amazon, to a continuous
#' application of CO~2~ of fossil origin. The main goal is to simulate a CO~2~ 
#' fertilization experiment analogous to Amazon FACE.
#' 
#' I use here the preliminary version of the carbon model developed for ATTO using
#' the `porce` modeling framework. Therefore, predictions include a reasonable level
#' of uncertainty, which is somehow consistent with the level of variation in carbon
#' stocks and radiocarbon fluxes measured at the terra firme forests in ATTO.
#' 
#' To model radiocarbon dynamics, it is necessary to model the incorporation of 
#' atmospheric radiocarbon in the ecosystem during the bomb period. Therefore, 
#' the simulations here are split in two steps. In the first step, I first model
#' the incorporation of bomb radiocarbon, from the year 1942 until 2024. In the
#' second step, I add radiocarbon with a signature of -1000 per mil in $\Delta^{14}$C
#' notation. This is equivalent to adding only radiocarbon free CO~2~ to the ecosystem. 
#' Simulations in this second step are done for two cases, a control case in which
#' only the label is added but there is no increase in GPP, and a second case in 
#' which elevated CO~2~ induce a 50% increase in GPP relative to the control case. 
#' These simulations are run for 10 years, assuming this as the time frame of the
#' experiment. 
#' 
#' ## Preparations
#' Let's first load all required packages, create some color palettes, fire up the
#' cluster and get a mean GPP for the site based on current observations. 
#' 
#+ message=FALSE
library(porce)
library(parallel)
library(SoilR)
library(RColorBrewer)
pal=brewer.pal(7, "Set2")
pal1<-rainbow(n=7, alpha=0.9)
pal2<-rainbow(n=7, alpha=0.5)
pal3<-hcl.colors(4, alpha=0.5)

no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)

setwd("~/MPIBGC-TEE/porce/simulation/ATTO/")
monthly_gpp<-read.csv("inputs/ATTO_INST_monthlyGPP_2014_2019.csv")
pool_names<-c("Foliage", "Wood", "Fine roots", "Coarse roots", "Fine litter", "CWD", "Soil (0-30 cm)")
GPPmean<- mean(unlist(monthly_gpp[-1])) * 1e-06 * 1e04 *12 
GPPsd<- sd(unlist(monthly_gpp[-1])) * 1e-06 * 1e04 *12 

set.seed(1033)
atto_gpp<-rnorm(n=nrow(modpars), mean=GPPmean, sd=GPPsd/sqrt(6))

#' We also need to create a time series of atmospheric radiocarbon for the bomb
#' period.

fNHZ3<-read.csv(file="inputs/NHZ3forecast.csv")
NHZ3<-data.frame(Year=c(Hua2021$NHZone3[-(937:939),1], fNHZ3$time),
                 Delta14C=c(Hua2021$NHZone3[-(937:939),2],Delta14C_from_AbsoluteFractionModern(fNHZ3$F14C)))

#' To solve the model over time, we need to specify initial conditions. We assume
#' that the forest is in dynamic equilibrium with respect to carbon. Therefore,
#' we create now a function that takes values of carbon stocks and some parameters
#' from the `porce` region/framework and calculates equilibrium carbon stocks.

pars2stocks<-function(gpp, pars){
  mod<-lam(input=inputGPP(gpp), matrix=makeB(pars))
  equilibriumStock(model=mod)
}

#' We use now this function to compute the initial C stocks for the entire set of 
#' parameter values (1000 sets). This is equivalent to obtaining 1000 initial 
#' conditions for the carbon stocks. This way we can quantify the uncertainty due
#' to parameter values in the simulations. 

allObjects=ls()
clusterExport(cl, varlist=c(allObjects, "lam", "inputGPP", "makeB", "equilibriumStock", "equilibriumOutflux"))

stockList<-clusterMap(cl, fun=pars2stocks, gpp=atto_gpp, pars=asplit(modpars, MARGIN = 1))
stocksATTO<-matrix(unlist(stockList), nrow=7, ncol=1000)

#' Time depend simulations are run in SoilR using the `Model_14` model class. 
#' In this case we are not interested in obtaining a time depend solution of the model,
#' but only the predicted values in 2025, which will be used as initial conditions
#' for the second step with the fertilization experiment. The function below
#' returns the radiocarbon content in all pools for a target year, in this case 
#' the year 2025.5. 

years<-seq(1942,2026, by=1/2) 
C14_ty<-function(gpp, pars, targetYear=2025.5){
  Mod<-Model_14(t=years,inputFluxes = inputGPP(gpp),A=makeB(pars), ivList=stocksATTO[,1],
                inputFc=BoundFc(NHZ3, lag=0, format="Delta14C"), 
                initialValF=ConstFc(c(-20, -25, -20, -25, -20, -30, -50), format="Delta14C"))
  C14<-getF14(Mod)
  C14[years==targetYear,]
}

#' The computation of these radiocarbon values consumes a significant among of
#' computing resources, so we run them in a cluster and later extract the mean
#' of the 1000 predictions. 

allObjects=ls()
clusterExport(cl, varlist=c(allObjects, "getF14", "ConstFc", "Model_14", "BoundFc", "inputGPP", "makeB"))
CS14<-clusterMap(cl,fun=C14_ty, gpp=atto_gpp, pars=asplit(x=modpars, MARGIN=1))
CS14_2025<-matrix(unlist(CS14), ncol=7, nrow=1000, byrow=TRUE)

mCS14<-apply(CS14_2025, MARGIN = 2, mean) # initial condition for experiment

#' ## Fertilization experiment
#' We will run simulations from 2025 until 2034 at monthly time steps. The function
#' below returns a constant value of radiocarbon of -1000 per mil for each
#' time step assuming that all radiocarbon in the air is of fossil origin. 

future<-seq(2025.5, 2035.5, by=1/12)
fAtmC14<-sapply(future, function(t) -1000 )

#' The following function computes the radiocarbon response of the ecosystem
#' for any value of gpp, model parameters, and initial C stocks. For the initial
#' radiocarbon values, we use those from the previous step. 

fossilResponse<-function(gpp, parms, iniStocks){
  mod<-Model_14(t=future,inputFluxes = inputGPP(gpp),A=makeB(parms), ivList=as.numeric(iniStocks),
                 inputFc=BoundFc(data.frame(future, fAtmC14), lag=0, format="Delta14C"), 
                 initialValF=ConstFc(mCS14, format="Delta14C"))
  C14pools<-getF14(mod)
  C14bulk<-getF14C(mod)
  C14resp<-getF14R(mod)
  Ct<-getC(mod)
  Rt<-getReleaseFlux(mod)
  
  return(list(C14p=C14pools, C14b=C14bulk, C14r=C14resp, Ct=Ct, Rt=Rt))
}

#' We run now all simulations in a cluster for all `modpars` values and initial
#' C stocks. Two simulations are prepared here, one with the average GPP, and 
#' another with 1.5 times higher GPP. The last part of the code, compiles the results
#' and aggregates them according to mean and standard deviation. 
#' 
allObjects=ls()
clusterExport(cl, varlist=c(allObjects, "getF14", "getF14C", "getF14R", "ConstFc", "Model_14", "BoundFc", "inputGPP", "makeB", "getC", "getReleaseFlux"))
ctrl<-clusterMap(cl, fun=fossilResponse, gpp=GPPmean, parms=asplit(modpars, 1), iniStocks = asplit(stocksATTO, 2))

ctrlC14b<-sapply(ctrl, FUN=function(x) x$C14b)
ctrlC14bsts<-data.frame(mean=apply(ctrlC14b, 1, mean), sd=apply(ctrlC14b, 1, sd))
ctrlC14r<-sapply(ctrl, FUN=function(x) x$C14r)
ctrlC14rsts<-data.frame(mean=apply(ctrlC14r, 1, mean), sd=apply(ctrlC14r, 1, sd))
ctrlC14p<-array(unlist(lapply(ctrl, FUN=function(x) x$C14p)), dim=c(length(future), 7, 1000))
ctrlC14psts<-list(mean=apply(ctrlC14p, c(1,2), mean), sd=apply(ctrlC14p, c(1,2), sd))
ctrlCt<-array(unlist(lapply(ctrl, FUN=function(x) x$Ct)), dim=c(length(future), 7, 1000))
ctrlCtsts<-list(mean=apply(ctrlCt, c(1,2), mean), sd=apply(ctrlCt, c(1,2), sd))
ctrlRt<-array(unlist(lapply(ctrl, FUN=function(x) x$Rt)), dim=c(length(future), 7, 1000))
ctrlRtsts<-list(mean=apply(ctrlRt, c(1,2), mean), sd=apply(ctrlRt, c(1,2), sd))


dbl<-clusterMap(cl, fun=fossilResponse, gpp=GPPmean*1.5, parms=asplit(modpars, 1), iniStocks = asplit(stocksATTO, 2))

dblC14b<-sapply(dbl, FUN=function(x) x$C14b)
dblC14bsts<-data.frame(mean=apply(dblC14b, 1, mean), sd=apply(dblC14b, 1, sd))
dblC14r<-sapply(dbl, FUN=function(x) x$C14r)
dblC14rsts<-data.frame(mean=apply(dblC14r, 1, mean), sd=apply(dblC14r, 1, sd))
dblC14p<-array(unlist(lapply(dbl, FUN=function(x) x$C14p)), dim=c(length(future), 7, 1000))
dblC14psts<-list(mean=apply(dblC14p, c(1,2), mean), sd=apply(dblC14p, c(1,2), sd))
dblCt<-array(unlist(lapply(dbl, FUN=function(x) x$Ct)), dim=c(length(future), 7, 1000))
dblCtsts<-list(mean=apply(dblCt, c(1,2), mean), sd=apply(dblCt, c(1,2), sd))
dblRt<-array(unlist(lapply(dbl, FUN=function(x) x$Rt)), dim=c(length(future), 7, 1000))
dblRtsts<-list(mean=apply(dblRt, c(1,2), mean), sd=apply(dblRt, c(1,2), sd))

pol<-function(x, df, cl){
  polygon(c(x, rev(x)), c(df$mean+df$sd, rev(df$mean - df$sd)),
          border=cl, col=cl)
}

#' ## Results
#' The results show that increasing GPP with a fossil radiocarbon signal leads to
#' a faster incorporation of fossil radiocarbon carbon in comparison with the 
#' control treatment. The response is faster in the ecosystem respiration flux
#' than in the ecosystem carbon pools. In fact, the respiration flux reaches a
#' quasi-equilibrium with the atmosphere. In contrast, after 10 years of fertilization
#' the ecosystem pools are still responding to the pulse and have not reached 
#' isotopic equilibration with the atmosphere. 

#+ fig.width=10, fig.height=8
pdf("~/SOIL-R/Meetings/2025/ATTOworkshop/Fert14C.pdf", encoding = "WinAnsi.enc")
par(mar=c(4,4.5,0,1))
plot(NULL, type="l", xlab="Year AD", ylab=expression(paste(Delta^14, "C (\u2030)")), 
     xlim=c(2025, 2035), ylim=c(-1000, 100), bty="n")
pol(future, ctrlC14bsts, pal3[1])
pol(future, ctrlC14rsts, pal3[2])
pol(future, dblC14bsts, pal3[3])
pol(future, dblC14rsts, pal3[4])
legend("bottomleft", c("Control, ecosystem radiocarbon", "Control, respired radiocarbon",
                       "+50%GPP, ecosystem radiocarbon", "+50%GPP, respired radiocarbon"), lty=1, lwd=4, col=pal3, bty="n")
dev.off()

#' The response of individual pools shows that fertilization with fossil carbon
#' leads to faster incorporation of radiocarbon in all pools (lighter color). However,
#' due to parameter uncertainty, this effect is easier to see in some pools than in others.
#' 

#+ fig.width=10, fig.height=12
pdf("~/SOIL-R/Meetings/2025/ATTOworkshop/Pools14C.pdf", height=7*sqrt(2), encoding = "WinAnsi.enc")
par(mfrow=c(4,2))
for(i in 1:7){
  plot(NULL, type="l", xlab="Year AD", ylab=expression(paste(Delta^14, "C (\u2030)")), 
       xlim=c(2025, 2035), ylim=c(-1000, 100), bty="n", main=pool_names[i])
    pol(future, df=data.frame(mean=ctrlC14psts$mean[,i], sd=ctrlC14psts$sd[,i]), cl=pal1[i])
    pol(future, df=data.frame(mean=dblC14psts$mean[,i], sd=dblC14psts$sd[,i]), cl=pal2[i])

}
par(mfrow=c(1,1))
dev.off()

#' The response of individual pools in terms of C stocks also suggest a large
#' degree of variability and responses due to fertilization are difficult to see,
#' except for the carbon stock in foliage and in the fine litter pool. 

#+ fig.width=10, fig.height=12
pdf("~/SOIL-R/Meetings/2025/ATTOworkshop/PoolsC.pdf", height=7*sqrt(2), encoding = "WinAnsi.enc")
par(mfrow=c(4,2))
for(i in 1:7){
  plot(NULL, type="l", xlab="Year AD", ylab=expression(paste("Carbon stock (Mg C h", a^-1, " y", r^-1, ")")), 
       xlim=c(2025, 2035), ylim=c(0, max(dblCtsts$mean[,i])*1.25), bty="n", main=pool_names[i])
  pol(future, df=data.frame(mean=ctrlCtsts$mean[,i], sd=ctrlCtsts$sd[,i]), cl=pal1[i])
  pol(future, df=data.frame(mean=dblCtsts$mean[,i], sd=dblCtsts$sd[,i]), cl=pal2[i])
  
}
par(mfrow=c(1,1))
dev.off()

#' The total response of the ecosystem in terms of the total stock overlaps
#' significantly with the background variability of the control treatment. 
#' If total stocks would be possible to quantify in the field, statistical 
#' differences among treatments may not be possible to see for the first 
#' 10 years of the experiment.

#+ fig.width=8, fig.height=10

pdf("~/SOIL-R/Meetings/2025/ATTOworkshop/TotalC.pdf", encoding = "WinAnsi.enc")
plot(NULL, type="l", xlab="Year AD", ylab=expression(paste("Total carbon stock (Mg C h", a^-1, ")")), 
     xlim=c(2025, 2035), ylim=c(200, 500), bty="n")
pol(future, df=data.frame(mean=rowSums(ctrlCtsts$mean), sd=sqrt(rowSums(ctrlCtsts$sd^2))), cl=pal3[1])
pol(future, df=data.frame(mean=rowSums(dblCtsts$mean), sd=sqrt(rowSums(dblCtsts$sd^2))), cl=pal3[2])
legend("topleft", c("Control", "+50% GPP"), pch=15, col=pal3, bty="n")
dev.off()

#' Aboveground biomass carbon is one of the most common measurements in tropical forests.
#' Model predictions also suggests that statistical differences may not be possible
#' to detect in this ecosystem component.

#+ fig.width=8, fig.height=10
abvg_pools<-match(c("Foliage", "Wood"), table=pool_names)
plot(NULL, type="l", xlab="Year AD", ylab=expression(paste("Aboveground carbon stock (Mg C h", a^-1, ")")), 
     xlim=c(2024, 2034), ylim=c(0, 200), bty="n")
pol(future, df=data.frame(mean=rowSums(ctrlCtsts$mean[,abvg_pools]), sd=sqrt(rowSums(ctrlCtsts$sd[,abvg_pools]^2))), cl=pal3[1])
pol(future, df=data.frame(mean=rowSums(dblCtsts$mean[,abvg_pools]), sd=sqrt(rowSums(dblCtsts$sd[,abvg_pools]^2))), cl=pal3[2])

#' The model does predict important differences in terms of ecosystem respiration.
#' This suggest that extra carbon entering through photosynthesis is mostly used
#' for respiration.

#+ fig.width=8, fig.height=10
pdf("~/SOIL-R/Meetings/2025/ATTOworkshop/ER.pdf")
plot(NULL, type="l", xlab="Year AD", ylab=expression(paste("Ecosystem respiration (Mg C h", a^-1, " y", r^-1, ")")), 
     xlim=c(2025, 2035), ylim=c(20, 60), bty="n")
pol(future, df=data.frame(mean=rowSums(ctrlRtsts$mean), sd=sqrt(rowSums(ctrlRtsts$sd^2))), cl=pal3[1])
pol(future, df=data.frame(mean=rowSums(dblRtsts$mean), sd=sqrt(rowSums(dblRtsts$sd^2))), cl=pal3[2])
legend("topleft", c("Control", "+50% GPP"), pch=15, col=pal3, bty="n")
dev.off()

#' In particular, autotrophic respiration is responsible for most of the response. 

#+ fig.width=8, fig.height=10
autotroph_pools=match(c("Foliage", "Wood", "Fine roots", "Coarse roots"), table=pool_names)
heterotroph_pools=match(c("Fine litter", "CWD", "Soil (0-30 cm)"), table=pool_names)

pdf("~/SOIL-R/Meetings/2025/ATTOworkshop/AR.pdf")
plot(NULL, type="l", xlab="Year AD", ylab=expression(paste("Autotrophic respiration (Mg C h", a^-1, " y", r^-1, ")")), 
     xlim=c(2025, 2035), ylim=c(10, 40), bty="n")
pol(future, df=data.frame(mean=rowSums(ctrlRtsts$mean[,autotroph_pools]), sd=sqrt(rowSums(ctrlRtsts$sd[,autotroph_pools]^2))), cl=pal3[1])
pol(future, df=data.frame(mean=rowSums(dblRtsts$mean[,autotroph_pools]), sd=sqrt(rowSums(dblRtsts$sd[,autotroph_pools]^2))), cl=pal3[2])
legend("topleft", c("Control", "+50% GPP"), pch=15, col=pal3, bty="n")
dev.off()

pdf("~/SOIL-R/Meetings/2025/ATTOworkshop/HR.pdf")
plot(NULL, type="l", xlab="Year AD", ylab=expression(paste("Heterotrophic respiration (Mg C h", a^-1, " y", r^-1, ")")), 
     xlim=c(2025, 2035), ylim=c(0, 20), bty="n")
pol(future, df=data.frame(mean=rowSums(ctrlRtsts$mean[,heterotroph_pools]), sd=sqrt(rowSums(ctrlRtsts$sd[,heterotroph_pools]^2))), cl=pal3[1])
pol(future, df=data.frame(mean=rowSums(dblRtsts$mean[,heterotroph_pools]), sd=sqrt(rowSums(dblRtsts$sd[,heterotroph_pools]^2))), cl=pal3[2])
legend("topleft", c("Control", "+50% GPP"), pch=15, col=pal3, bty="n")
dev.off()
