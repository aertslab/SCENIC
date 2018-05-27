# Saves the list of TFs and the transposed expression matrix as text, to load from GRNBoost
#' @title exportsForGRNBoost
#' @description Exports the TF list and expression matrix to text format to run GRNBoost (https://arboretum.readthedocs.io).
#' @param exprMat Expression matrix
#' @param scenicOptions Fields used: getDatasetInfo(scenicOptions, "org"), to select the annotation database.
#' @return Writes the output as text files.
#' @examples 
#' exportsForGRNBoost(exprMat, scenicOptions)
#' @export
exportsForGRNBoost <- function(exprMat, scenicOptions)
{
  allTFs <- getDbTfs(scenicOptions)
  
  inputTFs <- allTFs[allTFs %in% rownames(exprMat)] 
  write.table(inputTFs, file="int/1.2_inputTFs.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  message("TF list written as: ", "int/1.2_inputTFs.txt")
  
  
  exprMat_t <- t(exprMat)
  write.table(exprMat_t, file="int/1.1_exprMatrix_filtered_t.txt", sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)
  message("Transposed expression matrix written as: ", "int/1.1_exprMatrix_filtered_t.txt")
  
  
  # linkList <- GENIE3::getLinkList(weightMatrix, threshold=getSettings(scenicOptions, "modules/weightThreshold"))
  # colnames(linkList) <- c("TF", "Target", "weight")
  # 
  # # Order by weight
  # linkList <- linkList[order(linkList[,"weight"], decreasing=TRUE),]
  # saveRDS(linkList, file=getIntName(scenicOptions, "genie3ll"))
}


# TODO
#' @title importGRNBoostResults
#' @description Imports and formats GRNBoost results to continue the SCENIC workflow
#' @param scenicOptions Fields used: 
#' @return # TODO
#' @examples 
#' importGRNBoostResults(scenicOptions)
# export
importGRNBoostResults <- function(scenicOptions, ...)
{
 # options Wegight threshold <- XXX 
  
  
  # With GRNboost  # TODO
  # import
  # colnames (in order) c("TF", "Target", "weight")
  # save
  
  
  # linkList <- GENIE3::getLinkList(weightMatrix, threshold=getSettings(scenicOptions, "modules/weightThreshold"))
  # colnames(linkList) <- c("TF", "Target", "weight")
  # 
  # # Order by weight
  # linkList <- linkList[order(linkList[,"weight"], decreasing=TRUE),]
  # saveRDS(linkList, file=getIntName(scenicOptions, "genie3ll"))
}

