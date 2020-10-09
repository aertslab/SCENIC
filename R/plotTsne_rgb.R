#' @title plotTsne_rgb
#' @description Colors the t-SNE based on the activity of 3 (groups of) regulons
#' @param scenicOptions Fields used: AUC matrix (continuous or binary), default t-SNE.
#' @param regulonNames Regulons to plot
#' @param aucType "AUC" or "Binary"
#' @param aucMaxContrast To increase the AUC contrast decrease the value.
#' @param offColor Color por the cells completelly off. To deactivate (color as black), set to NULL.
#' @param showPlot Whether to plot the coloured t-SNE.
#' @param showPlot Whether to plot add a legend to the plot.
#' @param tSNE_fileName tSNE file name. If null, the default t-SNE is used.
#' @param ... Other arguments to pass to the \code{plot} function.
#' @return The cell colors (invisible)
#' @examples 
#' par(mfrow=c(1,2))
#' 
#' regulonNames <- c("Dlx1", "Sox8")
#' cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
#' text(-5,-23, attr(cellCol,"red"), col="red", cex=.7)
#' text(-10,-18, attr(cellCol,"green"), col="green", cex=.7)
#' 
#' regulonNames <- list(red=c("Dlx1","Dlx5"),
#'                      green=c("Neurod1"),
#'                      blue=c("Sox8"))
#' cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="Binary")
#'
#' @export
plotTsne_rgb <- function(scenicOptions, regulonNames, aucType="AUC", 
                         aucMaxContrast=0.8, offColor="#c0c0c030", 
                         showPlot=TRUE, showLegend=TRUE, tSNE_fileName=NULL, ...)
{
  # Check format
  if(length(regulonNames)>3) stop("To plot more than three regulons, group them by color.")
  if(is.null(names(regulonNames)) && length(regulonNames)<=3) names(regulonNames) <- c("red","green", "blue")[seq_along(regulonNames)]
  if(any(!names(regulonNames) %in% c("red","green", "blue"))) 
      stop('If a list, the names of regulonNames should be "red","green", and/or"blue".')
  if(aucMaxContrast<=0 | aucMaxContrast>1) stop("aucMaxContrast should be between 0 and 1")
    
  # Load data
  aucType <- tolower(aucType)
  if(aucType=="binary") mat4col <- loadInt(scenicOptions, "aucell_binary_nonDupl") 
  if(aucType=="auc") mat4col <- getAUC(loadInt(scenicOptions, "aucell_regulonAUC"))
  if(showPlot && is.null(tSNE_fileName)) tSNE_fileName <- tsneFileName(scenicOptions)
  tSNE <- readRDS(tSNE_fileName)
  
  mat4col <- mat4col[onlyNonDuplicatedExtended(rownames(mat4col)),,drop=FALSE]
  reguCols <- lapply(regulonNames, function(x) unlist(sapply(x, function(tf) {
    rownames(mat4col)[which(gsub("_extended","",getTF(rownames(mat4col))) %in% tf)] 
    }
    )))
  reguCols <- reguCols[lengths(reguCols)>0]
  
  # Average of binary...
  if(aucType=="binary")
    cellColChan <- sapply(reguCols, function(modsCol) apply(mat4col[modsCol,, drop=FALSE], 2, mean))
  # Color if all modules are "on"
  # cellColChan <- sapply(reguCols, function(modsCol) apply(mat4col[modsCol,, drop=FALSE], 2, function(x) as.numeric(sum(x)==length(x))))
  
  # AUC
  if(aucType=="auc") {
    cellColChan <- sapply(reguCols, function(regCol) {
      aucAvg <- colMeans(mat4col[regCol,, drop=FALSE])
      setNames(sapply(as.numeric(aucAvg/(max(aucAvg)*aucMaxContrast)), min, 1), names(aucAvg))
    })
  }
    
  
  # Apply color
  missingCol <- setdiff(c("red","green", "blue"), colnames(cellColChan))
  if(length(missingCol)>0)
    cellColChan <- cbind(cellColChan, matrix(rep(0, nrow(cellColChan)*length(missingCol)),ncol=length(missingCol), dimnames=list(NULL,missingCol)))
  cellCol <- apply(cellColChan, 1, function(x) rgb(x["red"], x["green"], x["blue"], alpha=.8))
  if(!is.null(offColor)) cellCol[which(cellCol=="#000000CC")] <- offColor # mostly for binary
  names(cellCol) <- colnames(mat4col)
  
  if(aucType=="binary") attr(cellCol,"Description") <- "Color: average of BINARY regulon activity"
  if(aucType=="auc") attr(cellCol,"Description") <- "Color: average of regulon activity (AUC)"
  attr(cellCol, "red") <- paste0(reguCols$red, collapse=", ")
  attr(cellCol, "green") <- paste0(reguCols$green, collapse=", ")
  attr(cellCol, "blue") <- paste0(reguCols$blue, collapse=", ")
  
  if(showPlot) 
  {
    plot(tSNE$Y, col=cellCol[rownames(tSNE$Y)], pch=16, 
         sub=attr(cellCol, "Description"), axes=FALSE, ...)
    
    if(showLegend)
    {
      cellColChan[which(cellColChan < (aucMaxContrast/2), arr.ind = T)] <- 0
      cellColChan <- cellColChan[which(apply(cellColChan, 1, function(x) any(x>0))),]
      cellLabels <- setNames(colnames(cellColChan)[apply(cellColChan, 1, function(x) which.max(x))], rownames(cellColChan))
      labsCoords <- t(sapply(split(data.frame(tSNE$Y), as.character(cellLabels[rownames(tSNE$Y)])), colMeans));
      
      regulonNames[rownames(labsCoords)] <- sapply(regulonNames[rownames(labsCoords)], paste, collapse=", ")
      for(i in rownames(labsCoords)) text(mean(labsCoords[i,1]), mean(labsCoords[i,2]), regulonNames[i], cex=1, col="black")
    }
  }
  invisible(cellCol)
}

