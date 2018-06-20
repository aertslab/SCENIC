#' @title export2scope
#' @description Create .loom file with the results of the analysis to visualize them in http://scope.aertslab.org
#' @param scenicOptions scenicOptions object
#' @param dgem Digital expression matrix
#' @return The .loom file (\code{file name indicated in "loomFile" slot in scenicOptions})containing the following information:
#' \itemize{
#' \item dgem (digital expression matrix)
#' \item title: \code{getDatasetInfo(scenicOptions,"datasetTitle")}
#' \item genome:\code{ getDatasetInfo(scenicOptions,"org")}
#' \item regulons.AUC: \code{getAUC(loadInt(scenicOptions, "aucell_regulonAUC"))}
#' \item regulons: \code{loadInt(scenicOptions, "aucell_regulons")}
#' \item default.embedding: default t-SNE
#' \item default.embedding.name: default t-SNE type (e.g. AUC/binary and description)
#' \item column attributes: \code{getDatasetInfo(scenicOptions, "cellInfo")} and nGene (\code{colSums(dgem>0)}).
#' }
#' @seealso To add extra data (e.g. embeddings or clusters), see: \code{help(package="SCopeLoomR")}
#' @examples 
#' # devtools::install_github("aertslab/SCopeLoomR")
#' # scenicOptions@fileNames$output["loomFile",] <- "myAnalysis.loom"
#' # export2scope(scenicOptions, exprMat)
#' @export 
export2scope <- function(scenicOptions, dgem)
{
  # TODO: what about if there is no dgem, but only normalized?
  # TODO: Error File 'int/cellInfo.Rds' does not exist.
  
  
  # TODO: ask max about order of samples tsne-expr-info
  suppressPackageStartupMessages(require(SCopeLoomR))

  # Default embedding (e.g. t-SNE or PCA coordinates)
  defaultTsne <- readRDS(tsneFileName(scenicOptions))
  defaultTsne_name <- paste("SCENIC t-SNE:", defaultTsne$type)
  defaultTsne <- defaultTsne$Y
  
  fileName <- getOutName(scenicOptions, "loomFile")
  if(file.exists(fileName)) stop("File '", getOutName(scenicOptions, "loomFile"), "' already exists.")
  loom <- build_loom(file.name=fileName,
                     dgem=dgem,
                     title=getDatasetInfo(scenicOptions,"datasetTitle"),
                     genome=getDatasetInfo(scenicOptions,"org"), # Just for user information, not used internally
                     default.embedding=defaultTsne,
                     default.embedding.name=defaultTsne_name)
  
  # Known cell information/annotation  
  cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
  if(!is.null(cellInfo))
  {
    cellInfo <- data.frame(cellInfo)
    cellInfo <- cellInfo[,colnames(cellInfo) != "nGene", drop=FALSE]
  }
  
  # Add annotation
  for(cn in colnames(cellInfo))
  {
    add_col_attr(loom=loom, key=cn, value=cellInfo[,cn], as.md.annotation=T) # as.md.annotation: to plot on tSNE
  }
  cellInfo$nGene <- colSums(dgem>0)
  add_col_attr(loom=loom, key = "nGene", value=cellInfo$nGene)
  
  
  # Regulons AUC matrix
  regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
  add_scenic_regulons_auc_matrix(loom=loom, regulons.AUC=getAUC(regulonAUC))
  
  # Regulons (gene list)
  regulons <- loadInt(scenicOptions, "aucell_regulons")
  add_scenic_regulons(loom=loom, dgem=dgem, regulons=regulons)
  
  # Add regulons thresholds #TODO
  # add_global_md_regulon_thresholds(loom=loom, regulon.threshold.assignments=...$threshold.assignment)
  
  # # Alternative t-SNE #TODO
  # load(paste0(scenicDir, "int/3.3_tsneRegulonAUC_PCA_50pc50perpl.RData"))
  # add_embedding(loom=loom, embedding=tSNE$Y, name="SCENIC (t-SNE on AUC)")
  
  finalize(loom=loom)
}
