#' @title plotTsne_AUCellApp
#' @description Creates the AUCell shiny app to visualize the regulon activities on the t-SNE
#' @param scenicOptions Fields used: aucell_regulonAUC, aucell_thresholds, cellInfo, and default tSNE
#' @param exprMat Expression matrix
#' @param tSNE_fileName tSNE file name. If null, the default t-SNE is used.
#' @param skipDuplicatedExtended Skip extended regulons if there is a high-confidence regulon for the same TF
#' @return Invisible: aucellApp 
#' @examples 
#' logMat <- exprMat # Better if it is logged/normalized
#' aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
#' savedSelections <- shiny::runApp(aucellApp)
#' @export
plotTsne_AUCellApp <- function(scenicOptions, exprMat, tSNE_fileName=NULL, skipDuplicatedExtended=TRUE) 
{
  if(is.null(tSNE_fileName)) tSNE_fileName <- tsneFileName(scenicOptions)
  tSNE <- readRDS(tSNE_fileName)
  regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
  thresholds <- loadInt(scenicOptions, "aucell_thresholds", ifNotExists="null")
  cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
  colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
  
  if(skipDuplicatedExtended)
  {
    regulons2plot <- onlyNonDuplicatedExtended(rownames(regulonAUC))
    regulonAUC <- regulonAUC[regulons2plot,]
    thresholds <- thresholds[regulons2plot]
  }
  aucellApp <- AUCell_createViewerApp(auc=regulonAUC,
                                      thresholds=thresholds,
                                      tSNE=tSNE$Y, 
                                      exprMat=exprMat,
                                      cellInfo=cellInfo,
                                      colVars=colVars)
  invisible(aucellApp)
}