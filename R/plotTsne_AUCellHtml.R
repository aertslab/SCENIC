################################################################
# Plot AUC histograms & default tsne

#' @title plotTsne_AUCellHtml
#' @description Wrapper for the 'AUCell_plotTSNE' function: plots the regulon activities on the t-SNE.
#' @param scenicOptions Fields used: aucell_regulonAUC, aucell_thresholds (int), and if \code{tSNE=NULL} the default tSNE.
#' @param exprMat Expression matrix (if NULL, TF expression will not be plotted)
#' @param fileName Location for the plots.
#' @param tSNE If null, it plots the default t-SNE.
#' @return NULL. The plots are saved as HTML with the name stored in \code{fileName}. 
#' If you are using RStudio, you might need to download the .html file and plots folder to view locally.
#' @examples 
#' plotTsne_AUCellHtml(scenicOptions, exprMat, fileName)
#' @export
plotTsne_AUCellHtml <- function(scenicOptions, exprMat, fileName, tSNE=NULL) 
{
  regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
  if(is.null(tSNE)) tSNE <- readRDS(tsneFileName(scenicOptions))
  tSNE <- tSNE$Y

  # Plot module activity, thresholds & assignment:
  thresholds <- loadInt(scenicOptions, "aucell_thresholds", ifNotExists="null")
  
  if(!is.null(thresholds)) 
    plots <- c("histogram", "binaryAUC", "AUC", "expression")
  if(is.null(thresholds)) 
    plots <- c("histogram", "AUC", "expression")

  plotsLoc <- dirname(fileName)
  plotsName <- basename(fileName)
  AUCell_plotTSNE(tSNE=tSNE, exprMat=exprMat,
                  cellsAUC=regulonAUC, thresholds=thresholds,
                  plots = c("histogram", "binaryAUC", "AUC", "expression"),
                  asPNG=plotsName) #alphaOff=0.1)
  file.rename(plotsName, file.path(plotsLoc, plotsName))

  if(file.exists("index.html"))
  {
    file.rename("index.html",
                file.path(plotsLoc, paste0(plotsName, ".html")))
  }else{
    if(!"R2HTML" %in% installed.packages()) warning("R2HTML package is not installed. Cannot produce HTML overview.")
  }
  
  invisible(NULL)
}