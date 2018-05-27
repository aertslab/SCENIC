# TODO: Test

#' @title regulon_plotExpression
#' @description Plots the expression of the genes in the regulon (as heatmap)
#' @param exprMat Expression matrix
#' @param regulonsSelected Regulons to plot
#' @param nCells Number of cells to plot (random subset of cells for faster plotting)
#' @param cellInfo Cell labels to plot in the heatmap (data.frame)
#' @param colVars Colors for the cell labels
#' @param color Expression color range
#' @param ... Other arguments to pass to \code{NMF::aheatmap}
#' @return Plots the heatmap
#' @examples 
#' library(SingleCellExperiment)
#' load("data/sceMouseBrain.RData")
#' exprMat <- counts(sceMouseBrain) # normalized preferred
#' dim(exprMat)
#' 
#' cellInfo <- data.frame(loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo")))
#' colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"))
#' 
#' exprMat_log <- log2(exprMat+1)
#' regulons <- loadInt(scenicOptions, "regulons")
#' regulonNames <- getRegulonName(c("Dlx1") , names(regulons))
#' regulon_plotExpression(exprMat_log, regulons[regulonNames], cellInfo=cellInfo, colVars=colVars) #, filename="regulonExpression.pdf")
#' @export
regulon_plotExpression <- function(exprMat, regulonsSelected, nCells=500, 
        cellInfo=NULL, colVars=NULL, color=c("black","goldenrod","yellow"), ...)
{
  cells2plot <- sample(colnames(exprMat), nCells)
  if(!is.null(cellInfo)) cellInfo <- cellInfo[cells2plot,, drop=FALSE]

  regulons_df <- reshape2::melt(regulonsSelected)
  colnames(regulons_df) <- c("gene","regulon")

  regulons_df <- regulons_df[which(regulons_df[,"gene"]%in% rownames(exprMat)),]
  
  if(length(regulonsSelected)>1) rownames(regulons_df) <- paste0(regulons_df[,"regulon"],"__", regulons_df[,"gene"])
  if(length(regulonsSelected)==1) rownames(regulons_df) <-regulons_df[,"gene"]

  mat2plot <- exprMat[regulons_df[,"gene"],cells2plot]
  rownames(mat2plot) <- rownames(regulons_df)
  NMF::aheatmap(mat2plot,
                color=color,
                annCol=cellInfo,
                annColor=colVars,
                annRow=regulons_df[,"regulon",drop=FALSE], ...) # geneExpression
}
