#TODO: Convert regulons to GeneSet or a class?


#' @title regulonsToGeneSet
#' @description Converts the regulons stored as incidence matrix to gene sets (e.g. list of genes)
#' @param incidMat Incidence matrix with TFs as rows and GENES as columns
#' @return List of genes per regulon
#' @examples 
#' regulonsToGeneSet(incidMat)
#' @export 
regulonsToGeneLists <- function(incidMat)
{
  regulons <- list()
  for(i in rownames(incidMat))
    regulons[[i]] <-  colnames(incidMat)[which(incidMat[i,]==1)]
  return(regulons)
}

#' @title getRegulonName
#' @description Gets the regulon name for a given TF (only returns the 'extended' regulons if no directly-annotated regulon is available)
#' @param TFs Transcription factor name
#' @param allRegulonNames List of regulon names (e.g. rownames(AUC...))
#' @param sep Character used to separate the TF name, if different from space.
#' @return Named character vector
#' @details ...
#' @seealso List of vignettes included in the package: \code{vignette(package="SCENIC")}
#' @examples
#' reguNames <- c("Dlx1 (103g)", "Dlx1_extended (190g)", "Olig2_extended (29g)", "Sox9 (17g)")
#' getRegulonName("Dlx1", reguNames)
#' getRegulonName("Olig2", reguNames)
#' @export
getRegulonName <- function(TFs, allRegulonNames, sep=" ")
{
  if(all(grepl(paste0(sep,"\\("), allRegulonNames))) allRegulonNames <- sapply(strsplit(allRegulonNames, paste0(sep,"\\(")), function(x) x[1])
  
  ret <- sapply(setNames(TFs,TFs), function(TF) allRegulonNames[grep(paste0(TF, "$"), allRegulonNames)])
  ret <- c(ret, setNames(sapply(TFs, function(TF) allRegulonNames[grep(paste0(TF, "_extended"), allRegulonNames)]), TFs))
  ret <- unlist(ret)
  
  ret <- unlist(onlyNonDuplicatedExtended(ret))
  
  ret
}

#' @title getTF
#' @description Returns the TF associated to a regulon name (e.g. removes the #genes sufix)
#' @param regulonName Character vector containing regulon names
#' @param sep Character used to separate the TF name (as regexpr), if different from space.
#' @return Character vector containing the TF (including the "_extended" sufix if appropiate)
#' @details ...
#' @seealso List of vignettes included in the package: \code{vignette(package="SCENIC")}
#' @examples
#' reguNames <- c("Dlx1 (103g)", "Dlx1_extended (190g)", "Olig2_extended (29g)", "Sox9 (17g)")
#' getTF(reguNames)
#' @export
getTF <- function(regulonName, sep="\\s")
{
  gsub(paste0(sep,"\\(\\d+g)"), "",regulonName)
}

#' @title onlyNonDuplicatedExtended
#' @description Returns the regulon names filtering-out the "extended" regulons if there is a regulon based on high-confidence annotations
#' @param regulonNames Character vector containing the regulon names (e.g. rownames(AUC_))
#' @return Character vector
#' @seealso List of vignettes included in the package: \code{vignette(package="SCENIC")}
#' @examples
#' reguNames <- c("Dlx1 (103g)", "Dlx1_extended (190g)", "Olig2_extended (29g)", "Sox9 (17g)")
#' onlyNonDuplicatedExtended(reguNames)
#' @export
onlyNonDuplicatedExtended <- function(regulonNames)
{
  regulonNames <- unname(regulonNames)
  tfs <- getTF(regulonNames)
  tfs <- gsub("_extended", "", tfs)
  splitRegulons <- split(regulonNames, tfs)[unique(tfs)]
  
  ret <- sapply(splitRegulons, function(x) 
  {
    split(x, grepl("_extended", x))[[1]] # False (direct) will be first
  })
  
  return(ret)
}

#' @title selectRegulons
#' @description Selects the regulons for the given TFs
#' @param regulons Regulons (list) or regulonAUC object
#' @param tfs TFs to select
#' @param onlyNonDuplicatedExtended Whether to filter with the function 'onlyNonDuplicatedExtended'
#' @return The selected regulons
#' @examples
#' selectRegulons(regulons, c("Dlx5", "Olig2"))
#' selectRegulons(regulonAUC, "Tef")
#' @export
selectRegulons <- function(regulons, tfs, onlyNonDuplicatedExtended=FALSE)
{
  regulonNames <- names(regulons)
  
  tfs <- c(tfs, paste0(tfs, "_extended"))
  regulonNames <- setNames(getTF(regulonNames), regulonNames)
  selectedRegulons <- names(regulonNames[regulonNames %in%  tfs])
  if(onlyNonDuplicatedExtended) selectedRegulons <- unname(onlyNonDuplicatedExtended(selectedRegulons))
  
  if(class(regulons)=="aucellResults") ret <- regulons[selectedRegulons,]
  if(class(regulons)=="list") ret <- regulons[selectedRegulons]
  
  return(ret)
}

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
#' loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom")
#' exprMat <- SCopeLoomR::get_dgem(SCopeLoomR::open_loom(loomPath)) # normalized preferred
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
  cells2plot <- colnames(exprMat)
  if(nCells < ncol(exprMat)) cells2plot <- sample(cells2plot, nCells)
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

