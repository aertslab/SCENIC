
#' @title runGenie3
#' @description Runs GENIE3 and formats the output
#' @param exprMat Expression matrix
#' @param scenicOptions Fields used: genie3wm, genie3ll (int)
#' @param ... Other arguments to pass to GENIE3 (e.g. ntrees=500)
#' @return Writes the output as \code{getIntName(scenicOptions, "genie3ll")}
#' @examples 
#' library(SCENIC)
#' scenicOptions <- readRDS("int/scenicOptions.Rds")
#' 
#' load("data/sceMouseBrain.RData")
#' exprMat <- counts(sceMouseBrain)
#' 
#' genesKept <- loadInt(scenicOptions, "genesKept")
#' exprMatrix_filtered <- exprMat[genesKept,]
#' exprMat_filtered <- log2(exprMatrix_filtered+1) 
#' 
#' runGenie3(exprMat_filtered, scenicOptions)
#' @export 
runGenie3 <- function(exprMat, scenicOptions, nParts=10, ...)
{
  nCores <- getSettings(scenicOptions, "nCores")
  
  # Check expression matrix (e.g. not data.frame)
  if(is.data.frame(exprMat)) 
  {
    supportedClasses <- paste(gsub("AUCell_buildRankings,", "", methods("AUCell_buildRankings")), collapse=", ")
    supportedClasses <- gsub("-method", "", supportedClasses)
    
    stop("'exprMat' should be one of the following classes: ", supportedClasses, 
         "(data.frames are not supported. Please, convert the expression matrix to one of these classes.)")
  }
  if(any(table(rownames(exprMat))>1))
    stop("The rownames (gene id/name) in the expression matrix should be unique.s")
  
  # Check TFs
  allTFs <- getDbTfs(scenicOptions)
  
  inputTFs <- allTFs[allTFs %in% rownames(exprMat)] # TODO check
  percMatched <- length(inputTFs)/length(allTFs)
  if(getSettings(scenicOptions, "verbose")) message("Using ", length(inputTFs), " TFs as potential regulators...")
  if(percMatched < .40) warning("Only ", round(percMatched*100), "% of the ", length(allTFs)," TFs in the database were found in the dataset. Do they use the same gene IDs?")
  
  # Run on subsets of genes 
  # (dividing the original gene list into 10 pieces)
  genesSplit <- suppressWarnings(split(sort(rownames(exprMat)), 1:nParts))
  weightMatrices <- list()
  for(i in 1:length(genesSplit))
  {
    if(getSettings(scenicOptions, "verbose")) message("Running GENIE3 part ", i)
    set.seed(getSettings(scenicOptions, "seed"))
    weightMatrix <- GENIE3::GENIE3(exprMat, regulators=inputTFs, nCores=nCores, targets=genesSplit[[i]], ...)
    fileName <- gsub(".Rds$", paste0("_part_", i,".Rds"), getIntName(scenicOptions, "genie3wm"))
    saveRDS(weightMatrix, file=fileName)
    weightMatrices[[i]] <- weightMatrix
  }
  
  # Convert to Link list:
  linkList_list <- list()
  for(i in 1:length(genesSplit))
  {
    # If interrupted or run by pieces:
    # fileName <- gsub(".Rds$", paste0("_part_", i,".Rds"), getIntName(scenicOptions, "genie3wm"))
    # weightMatrix <- readRDS(fileName)
    weightMatrix <- weightMatrices[[i]]
    linkList_list[[i]] <- GENIE3::getLinkList(weightMatrix, threshold=getSettings(scenicOptions, "modules/weightThreshold"))
  }
  rm(weightMatrices)
  
  # Merge
  linkList <- do.call(rbind, linkList_list)
  colnames(linkList) <- c("TF", "Target", "weight")
  
  # Order by weight
  linkList <- linkList[order(linkList[,"weight"], decreasing=TRUE),]
  saveRDS(linkList, file=getIntName(scenicOptions, "genie3ll"))
  
  invisible(linkList)
}


