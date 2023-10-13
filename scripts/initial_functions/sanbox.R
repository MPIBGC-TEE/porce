
source("~/MPIBGC-TEE/porce/pkg/R/lam.R")
source("~/MPIBGC-TEE/porce/pkg/R/equilibriumStock.R")
source("~/MPIBGC-TEE/porce/pkg/R/equilibriumOutflux.R")

toyModel<-new("lam", input=c(1,2,3), matrix=diag(-1,3,3))
toyModel<-lam(input=c(1,2,3), matrix=diag(-1,3,3))

equilibriumStock(toyModel)
equilibriumOutflux(toyModel)

mlist<-list(toyModel, toyModel)
simplify2array(lapply(mlist, FUN=equilibriumStock))
