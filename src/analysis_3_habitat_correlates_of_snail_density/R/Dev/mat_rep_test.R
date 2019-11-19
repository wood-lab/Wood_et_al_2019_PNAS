

m1 <- matrix(1:4, nrow = 2, byrow = T)
x1 <- cbind(rnorm(100,0,1), rnorm(100,5,1))
x2 <- cbind(rnorm(100,10,1), rnorm(100,15,1))
y <- rnorm(100, x1[,1] * m1[1,1] + x1[,2] * m1[1,2] + x2[,1] * m1[2,1] + x2[,2] * m1[2,2], 1)



setwd("src")
library(TMB)
version = "mat_rep_test"
TMB::compile( paste0(version,".cpp") )
dyn.load(dynlib(version))
setwd("../")

Obj = TMB::MakeADFun( data = list( x1 = x1, x2 = x2, y = y ), parameters = list( m1 = m1, log_sigma = 0), checkParameterOrder=FALSE, DLL = version)
opt <- Optimize(Obj)
rep <- sdreport(Obj)
Obj$env$parameters
