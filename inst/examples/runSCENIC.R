# For a detailed example see the package vignette:
vignette("SCENIC_runningStep1nWrapper")


\dontrun{

setwd("SCENIC_MouseBrain")
load("data/esetMouseBrain.RData")
exprMat <- exprs(esetMouseBrain)
# Optional: add log for TF expression plot, it does not affect any other calculation
exprMat <- log2(exprMat+1)
dim(exprMat)

load("data/colVars.RData")
cellInfo <- pData(esetMouseBrain)[colnames(exprMat), names(colVars), drop=F]

library(SCENIC)
runSCENIC(exprMat=exprMat, org="mm9", cellInfo=cellInfo, colVars=colVars, nCores=4,  stepsToRun=c("1.2", "2", "3.1", "3.2", "4"))

}
