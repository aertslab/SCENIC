#' @title plotTsne_compareSettings
#' @import AUCell
#' @description Plots several t-SNEs coloured by cell properties (to compare different t-SNE settings)
#' @param fileNames t-SNE files to plot
#' @param scenicOptions imports cellInfo and colVars
#' @param varName If a varName is given, it plots all the t-SNEs (not saved to file), colouring only for that variable.
#' If null, it plots all the cellInfo colums for all the t-SNEs (each t-SNE in an individual file). 
#' To plot all cellInfo colums, but as individual plots, use \code{TRUE}.
#' @param ... Other arguments to pass to \code{AUCell::plotTsne_cellProps}
#' @return Plots
#' @examples
#' fileNames <- paste0("int/", grep("tSNE_", list.files("int"), value=T))
#' par(mfrow=c(2,2))
#' plotTsne_compareSettings(fileNames, scenicOptions, varName="CellType")
#' @export
plotTsne_compareSettings <- function(fileNames, scenicOptions, varName=NULL, ...)
{
  cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
  colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
  
  if(is.null(varName))
  {
    if(getSettings(scenicOptions, "verbose")) message("t-SNE plots saved as: ")
    for(tSNE_fileName in fileNames)
    {
      tSNE <- readRDS(tSNE_fileName)
      plot_fileName <- paste0(gsub(".Rds", "", tSNE_fileName), ".pdf") #in case it is different extension...
      
      sub <- ""; if("type" %in% names(tSNE)) sub <- paste0("t-SNE on ", tSNE$type)
      
      pdf(plot_fileName)
      AUCell::plotTsne_cellProps(tSNE$Y, cellInfo=cellInfo, colVars=colVars, sub=sub, ...)
      dev.off()
      if(getSettings(scenicOptions, "verbose")) message("  ", plot_fileName)
    }
  }
  if(!is.null(varName))
  {
    if(is.character(varName) & any(! varName %in% colnames(cellInfo))) 
      stop("Cannot find cell information named '", varName[which(!varName %in% colnames(cellInfo))]," (it should be a column of the 'cellInfo').")
    
    cellInfo <- cellInfo[,varName,drop=FALSE]
      if(ncol(cellInfo)==0) stop("No properties selected for plotting.")
    
    for(tSNE_fileName in fileNames)
    {
      tSNE <- readRDS(tSNE_fileName)
      sub <- ""; if("type" %in% names(tSNE)) sub <- paste0("t-SNE on ", tSNE$type)
      
      AUCell::plotTsne_cellProps(tSNE$Y, cellInfo=cellInfo, colVars=colVars, sub=sub, ...)
    }
  }
}
