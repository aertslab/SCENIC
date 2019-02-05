
#' @title geneFiltering
#' @description Gene filtering
#' @param exprMat Expression matrix
#' @param scenicOptions scenicOptions object (slots used: RcisTarget databases and genesKept)
#' @param minCountsPerGene Minimum counts per gene required
#' @param minSamples Minimum number of samples (cells) in which the gene should be detected
#' @examples 
#' setwd("SCENIC_MouseBrain")
#' scenicOptions <- readRDS("int/scenicOptions.Rds")
#' 
#' loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom")
#' exprMat <- SCopeLoomR::get_dgem(SCopeLoomR::open_loom(loomPath))
#' 
#' genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
#'                            minCountsPerGene=3*.01*ncol(exprMat),
#'                            minSamples=ncol(exprMat)*.01)
#' overrideOptions <- list(dbFilePath="databases/mm9-500bp-upstream-7species.mc9nr.feather", outFile_genesKept=NULL)
#' @export 
geneFiltering <- function(exprMat, scenicOptions,
                          minCountsPerGene=3*.01*ncol(exprMat),
                          minSamples=ncol(exprMat)*.01)
{
  # Load options: outFile_genesKept and dbFilePath
  outFile_genesKept <- NULL
  dbFilePath <- NULL
  if(class(scenicOptions) == "ScenicOptions")
  {
    dbFilePath <- getDatabases(scenicOptions)[[1]]
    outFile_genesKept <- getIntName(scenicOptions, "genesKept")
  }else{
    dbFilePath <- scenicOptions[["dbFilePath"]]
    outFile_genesKept <- scenicOptions[["outFile_genesKept"]]
  }
  if(is.null(dbFilePath)) stop("dbFilePath")
  
  # Check expression matrix (e.g. not factor)
  if(is.data.frame(exprMat)) 
  {
    supportedClasses <- paste(gsub("AUCell_buildRankings,", "", methods("AUCell_buildRankings")), collapse=", ")
    supportedClasses <- gsub("-method", "", supportedClasses)
    
    stop("'exprMat' should be one of the following classes: ", supportedClasses, 
         "(data.frames are not supported. Please, convert the expression matrix to one of these classes.)")
  }
  if(any(table(rownames(exprMat))>1))
    stop("The rownames (gene id/name) in the expression matrix should be unique.")
  
  # Calculate stats
  nCountsPerGene <- rowSums(exprMat, na.rm = T)
  nCellsPerGene <- rowSums(exprMat>0, na.rm = T)
  
  ## Show info
  message("Maximum value in the expression matrix: ", max(exprMat, na.rm=T))
  message("Ratio of detected vs non-detected: ", signif(sum(exprMat>0, na.rm=T) / sum(exprMat==0, na.rm=T), 2))
  message("Number of counts (in the dataset units) per gene:")
  print(summary(nCountsPerGene))
  message("Number of cells in which each gene is detected:")
  print(summary(nCellsPerGene))

  ## Filter
  message("\nNumber of genes left after applying the following filters (sequential):")
  # First filter
  # minCountsPerGene <- 3*.01*ncol(exprMat)
  genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minCountsPerGene)]
  message("\t", length(genesLeft_minReads), "\tgenes with counts per gene > ", minCountsPerGene)
  
  # Second filter
  # minSamples <- ncol(exprMat)*.01
  nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
  genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
  message("\t", length(genesLeft_minCells), "\tgenes detected in more than ",minSamples," cells")
  
  # Exclude genes missing from database:
  library(RcisTarget)
  motifRankings <- importRankings(dbFilePath) # either one, they should have the same genes
  genesInDatabase <- colnames(getRanking(motifRankings))
  
  genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
  message("\t", length(genesLeft_minCells_inDatabases), "\tgenes available in RcisTarget database")
  
  genesKept <- genesLeft_minCells_inDatabases
  if(!is.null(outFile_genesKept)){ 
    saveRDS(genesKept, file=outFile_genesKept)
    if(getSettings(scenicOptions, "verbose")) message("Gene list saved in ", outFile_genesKept)
  }
  return(genesKept)
}
