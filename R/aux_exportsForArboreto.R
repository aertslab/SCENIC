#' @title exportsForArboreto (alias: exportsForGRNBoost)
#' @description Exports the TF list and expression matrix to text format to run GRNBoost (https://arboreto.readthedocs.io).
#' Note that it is also possible to run the full pipeline of SCENIC in python. See https://pyscenic.readthedocs.io
#' @param exprMat Expression matrix
#' @param scenicOptions Fields used: getDatasetInfo(scenicOptions, "org"), to select the annotation database.
#' @param dir Output folder (by default 'int' in current working directory)
#' @return Writes the output as text files.
#' @examples 
#' exportsForArboreto(exprMat, scenicOptions)
#' @export

exportsForArboreto <- exportsForGRNBoost <- function(exprMat, scenicOptions, dir="int")
{
  allTFs <- getDbTfs(scenicOptions)
  
  inputTFs <- allTFs[allTFs %in% rownames(exprMat)] 
  write.table(inputTFs, file=file.path(dir,"1.1_inputTFs.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  message("TF list written as: ", file.path(dir,"1.1_inputTFs.txt"))
  
  exprMat_t <- t(exprMat)
  write.table(exprMat_t, file=file.path(dir,"1.1_exprMatrix_filtered_t.txt"), sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)
  message("Transposed expression matrix written as: ", file.path(dir,"1.1_exprMatrix_filtered_t.txt"))
}
