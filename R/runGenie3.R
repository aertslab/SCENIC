
#' @title runGenie3
#' @description Runs GENIE3 and formats the output
#' @param exprMat Expression matrix
#' @param scenicOptions Fields used: genie3wm, genie3ll (int)
#' @param nParts Number of pieces to fragment the partial results into (to check progress or in case the job gets interrupted). Default: 10
#' @param resumePreviousRun Reload partial results from a previous (interrupted) run and resume execution (default: FALSE). 
#' @param ... Other arguments to pass to GENIE3 (e.g. ntrees=500)
#' @return Writes the output as \code{getIntName(scenicOptions, "genie3ll")}
#' @examples 
#' library(SCENIC)
#' scenicOptions <- readRDS("int/scenicOptions.Rds")
#' 
#' loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom")
#' exprMat <- SCopeLoomR::get_dgem(SCopeLoomR::open_loom(loomPath))
#' 
#' genesKept <- loadInt(scenicOptions, "genesKept")
#' exprMatrix_filtered <- exprMat[genesKept,]
#' exprMat_filtered <- log2(exprMatrix_filtered+1) 
#' 
#' runGenie3(exprMat_filtered, scenicOptions)
#' @export 
runGenie3 <- function(exprMat, scenicOptions, nParts=10, resumePreviousRun=FALSE, ...)
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
    stop("The rownames (gene id/name) in the expression matrix should be unique.")
  
  # Check TFs
  allTFs <- getDbTfs(scenicOptions)
  
  inputTFs <- allTFs[allTFs %in% rownames(exprMat)] 
  percMatched <- length(inputTFs)/length(allTFs)
  if(getSettings(scenicOptions, "verbose")) message("Using ", length(inputTFs), " TFs as potential regulators...")
  if(percMatched < .40) warning("Only ", length(inputTFs) ," (", round(percMatched*100), "%) of the ", length(allTFs)," TFs in the database were found in the dataset. Do they use the same gene IDs?\n")

  # Run (on subsets of genes) ----
  # (dividing the original gene list into 10 pieces)
  weightMatrices <- list()
  ## If interrupted or run by pieces: reload partial results
  if(resumePreviousRun)
  {
    if(file.exists(getIntName(scenicOptions, "genie3ll")))
    {
      stop("The previous run already finished running (the output already exists: '",getIntName(scenicOptions, "genie3ll"),"'). \nTo re-run GENIE3 set 'resumePreviousRun=FALSE'.")
    }else{
      fileNames <- gsub(".Rds$", "_part_", getIntName(scenicOptions, "genie3wm"))
      intDir <- dirname(fileNames)
      fileNames <- list.files(dirname(fileNames), pattern=basename(fileNames))
      if(length(fileNames) > 0){
        for(fileName in fileNames)
        {
          i <- gsub(".Rds$", "", strsplit(fileName, "_part_")[[1]][[2]])
          weightMatrices[[i]] <- readRDS(file.path(intDir, fileName))
          # if(!all(genesSplit[[i]] %in% colnames(weightMatrices[[i]]))) 
          #   warning("The split of genes from the current expression matrix does not match the partial results loaded.") # should match, but just in case, lets re-run it
        }
      }
    }
  }
  genesDone <- unname(unlist(lapply(weightMatrices, colnames)))
  if(any(table(genesDone) > 1)) stop("Some genes are in several of the partial runs.")
  
  ## Run the remaining: 
  genesLeft <- setdiff(sort(rownames(exprMat)), genesDone)
  if(length(genesLeft) > 0)
  {
    partNames <- as.character(seq_len(nParts))
    partNames <- setdiff(partNames, names(weightMatrices))
    if(length(partNames)==0) {
      warning("Splitting the ", length(genesLeft), " genes left into ", nParts, " parts." )
      partNames <- as.character(max(as.numeric(names(weightMatrices))) + seq_len(nParts))
    }
    genesSplit <- suppressWarnings(split(genesLeft, partNames))
    for(i in names(genesSplit))
    {
      if(getSettings(scenicOptions, "verbose")) message("Running GENIE3 part ", i)
      set.seed(getSettings(scenicOptions, "seed"))
      weightMatrix <- GENIE3::GENIE3(exprMat, regulators=inputTFs, nCores=nCores, targets=genesSplit[[i]], ...)
      fileName <- gsub(".Rds$", paste0("_part_", i,".Rds"), getIntName(scenicOptions, "genie3wm"))
      saveRDS(weightMatrix, file=fileName)
      weightMatrices[[i]] <- weightMatrix
    }
  }
  
  # Convert to Link list:
  linkList_list <- list()
  for(i in names(weightMatrices))
  {
    weightMatrix <- weightMatrices[[i]]
    linkList_list[[i]] <- GENIE3::getLinkList(weightMatrix, threshold=getSettings(scenicOptions, "modules/weightThreshold"))
  }
  rm(weightMatrices)
  
  # Merge
  linkList <- do.call(rbind, linkList_list)
  colnames(linkList) <- c("TF", "Target", "weight")
  
  # Order by weight
  linkList <- linkList[order(linkList[,"weight"], decreasing=TRUE),]
  linkList <- unique(linkList)
  rownames(linkList) <- NULL 
  saveRDS(linkList, file=getIntName(scenicOptions, "genie3ll"))
  
  if(getSettings(scenicOptions, "verbose")) message("Finished running GENIE3.")
  invisible(linkList)
}


