

#' @title viewMotifs
#' @description Shows the motif information table as HTML using DT::datatable()
#' @param tableSubset Motif information table to show
#' @param motifCol Column name containing the motif logo ID (to pass to RcisTarget::addLogo). If NULL the logo is not added.
#' @param dbVersion Database version (to pass to RcisTarget::addLogo).
#' @param nSignif Number of significant digits to show for numeric columns
#' @param colsToShow (No warning is shown for missing columns)
#' @param options argument to pass to DT::datatable()
#' @param ... Any other other arguments to pass to DT::datatable()
#' @seealso \code{?DT::datatable()}
#' @return Writes the output in the file name stored in: \code{getIntName(scenicOptions, "corrMat")}
#' @examples 
#' library(SCENIC)
#' scenicOptions <- readRDS("int/scenicOptions.Rds")
#' 
#' regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
#' tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
#' viewMotifs(tableSubset)
#' 
#' motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
#' tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox8"]
#' viewMotifs(tableSubset)
#' @export 
viewMotifs <- function(tableSubset, 
                      motifCol=c("motif", "bestMotif", "MotifID"), dbVersion="v9",
                      nSignif=3,
                      colsToShow=c(motifEnrichment=c("motifDb", "logo", "NES", "geneSet", "TF_highConf"),
                                       regulonTargets=c("TF", "gene", "nMotifs", "bestMotif", "logo", "NES", "highConfAnnot", "Genie3Weight")),
                      options=list(pageLength=50), 
                      ...)
{
  # Add motif logo
  if(!is.null(motifCol)){
    motifCol <- motifCol[which(motifCol %in% colnames(tableSubset))]
    if(length(motifCol) == 1){
      tableSubset <- RcisTarget::addLogo(tableSubset, motifCol=motifCol, dbVersion=dbVersion, addHTML=TRUE)
      if(!is.null(colsToShow)) colsToShow <- c("logo", colsToShow)
    }else{
      stop("Please indicate the column containing the motif id (argument 'motifCol') or set it to NULL.")
    }
    # TODO 
    tableSubset <- tableSubset[grep("transfac_pro__", tableSubset[[motifCol]], invert = T),]
  }
 
  # For numeric columns, show only the number of significant digits...
  for(i in which(sapply(tableSubset, is.numeric))){
    tableSubset[[i]] <- signif(tableSubset[[i]], nSignif)
  }

  # Keep only requested columns
  if(!is.null(colsToShow)) {
    colsToShow <- unique(unname(unlist(colsToShow)))
    colsToShow <- colsToShow[which(colsToShow %in% colnames(tableSubset))]
    tableSubset <- tableSubset[, colsToShow, with=F]
  }

  # create...
  DT::datatable(tableSubset, 
              escape=FALSE, filter="top", options=options)#, ...)
}

#' @export
viewTable <- function(tableSubset, motifCol=NULL, nSignif=3, colsToShow=NULL, options=list(pageLength=50), ...)
{
  viewMotifs( tableSubset=tableSubset, motifCol=motifCol, nSignif=nSignif, colsToShow=colsToShow, options=options, ...)  
}
