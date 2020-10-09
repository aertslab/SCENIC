# 
# weightMatrices <- list()
# for(i in 1:10)
# {
#   weightMatrix <- readRDS(paste0("/ddn1/vol1/staging/leuven/stg_00002/lcb/saibar/Projects/packages/SCENIC_MouseBrain_2018/int/1.3_GENIE3_weightMatrix_part_",i,".Rds"))
#   weightMatrices[[i]] <- weightMatrix
# }
# weightMatrix <- do.call(cbind, weightMatrices)
# rm(weightMatrices)
# dim(weightMatrix)
# library(Matrix)
# # threshold=getSettings(scenicOptions, "modules/weightThreshold")
# threshold=0.001
# weightMatrix[which(weightMatrix<threshold, arr.ind=T)] <- 0
# sum(weightMatrix==0)/sum(weightMatrix>0)
# weightMatrix <- as(weightMatrix, "dgCMatrix")
# saveRDS(weightMatrix, file="weightMatrix.Rds")
# 
# library(GENIE3)
# ll <- getLinkList(as.matrix(weightMatrix), threshold = 0.001) # thr: otherwise keeps the zeroes
# nrow(ll)
# # Yes, it is the same:
# # gn3 <- readRDS("/ddn1/vol1/staging/leuven/stg_00002/lcb/saibar/Projects/packages/SCENIC_MouseBrain_2018/int/1.4_GENIE3_linkList.Rds"); nrow(gn3)
