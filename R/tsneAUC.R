#' @title tsneAUC
#' @description Calculates the t-SNE based on the regulon activity. (i.e. binary/continuous AUC)
#' @param scenicOptions Fields used: 
#' @param aucType "AUC" for continuous AUC, or "Binary" . By default it uses the value in \code{getSettings(scenicOptions,"defaultTsne/dims")}
#' @param nPcs Number of PCs to use for the t-SNE (can be several values, see examples). By default it uses the value in \code{getSettings(scenicOptions,"defaultTsne/dims")}
#' @param perpl Perplexity for the t-SNE (can be several values, see examples). By default it uses the value in \code{getSettings(scenicOptions,"aucell/tsne_perpl")}
#' @param filePrefix Prefix for the file to save the t-SNE as. It is saved with the following format: "PREFIX_aucType_nPcs_perpl.Rds"
#' @param seed Seed for the t-SNEs. By default it uses the value in \code{getSettings(scenicOptions,"seed")}
#' @param onlyHighConf By default (TRUE) it only uses only regulons based on high-confidence annotations (i.e. "extended" are skipped)
#' @param ... Other arguments to pass to the \code{Rtsne::Rtsne} function.
#' @return The file name(s) in whihc the tSNE is saved. 
#' @examples 
#' ## TO DO (See the vignette)
#' # tSNE <- readRDS(tsneFileName(scenicOptions))
#' # tsneAUC(scenicOptions)
#' @export
tsneAUC <- function(scenicOptions, aucType=NULL, nPcs=NULL, perpl=NULL, filePrefix=NULL, seed=NULL, onlyHighConf=FALSE, ...)
{
  if(is.null(filePrefix)) filePrefix <- getSettings(scenicOptions,"tSNE_filePrefix")
  if(is.null(aucType)) aucType <- getSettings(scenicOptions,"defaultTsne/aucType")
  if(is.null(nPcs)) nPcs <- getSettings(scenicOptions,"defaultTsne/dims")
  if(is.null(perpl)) perpl <- getSettings(scenicOptions,"defaultTsne/perpl")
  
  if(is.null(seed)) seed <- getSettings(scenicOptions,"seed")
  
  if(is.character(nPcs)) {
    if(any(is.na(as.integer(nPcs[which(tolower(nPcs)!="dist")])))) stop("nPcs should be integer or 'dist'" )
  }
  if(length(aucType)>1) stop('aucType should be either "Binary" or "AUC"')
  if(length(aucType)>1) stop('aucType should be either "Binary" or "AUC"')
  if(!aucType %in% c("Binary", "AUC")) stop('aucType should be "Binary" or "AUC"')
    
  if(aucType=="Binary") mat4tsne <- loadInt(scenicOptions, "aucell_binary_nonDupl") 
  if(aucType=="AUC") mat4tsne <- getAUC(loadInt(scenicOptions, "aucell_regulonAUC"))

  if(onlyHighConf)
  {
    mat4tsne_subset <- mat4tsne[grep("_extended",rownames(mat4tsne),invert=T, value=T),] 
  }else{
    mat4tsne_subset <- mat4tsne[onlyNonDuplicatedExtended(rownames(mat4tsne)),]
  }
  
  # Prepare all possible combinations...
  allParams <- list()
  for(perplexity in perpl)
  {
    for(initial_dims in nPcs)
    {
      allParams[[length(allParams)+1]] <- c(perplexity=perplexity, initial_dims=initial_dims)
    }
  }
  
  if(getSettings(scenicOptions, "nCores")>1){
    suppressMessages(require("doMC", quietly=TRUE))
    doMC::registerDoMC(cores=getSettings(scenicOptions, "nCores"))
  }
 
  fileNames <- c()
  fileNames <- foreach(param=allParams) %dopar%
  {
    initial_dims <- param["initial_dims"]
    perplexity <- param["perplexity"]
    dimsAsText <- initial_dims
    if(tolower(initial_dims)!="dist") dimsAsText <- paste0(dimsAsText, "PCs")
    #if(getSettings(scenicOptions, "verbose")) message(paste0(format(Sys.time(), "%H:%M"), "\tCalculating t-SNE (", perplexity,"perpl, ", tolower(dimsAsText),")"))
    
    tryCatch(
    {
      # PCA-based t-SNE
      set.seed(seed)
      if(tolower(initial_dims)!="dist")
      {
        tsneMat <- mat4tsne_subset
        if(aucType=="Binary") tsneMat <- jitter(tsneMat, factor=1)
        tsneMat <- unique(t(tsneMat))
  
        tsneAUC <- Rtsne::Rtsne(tsneMat,
                                initial_dims=as.integer(initial_dims),
                                perplexity=perplexity, ...) #TODO tSNE function?
        rownames(tsneAUC$Y) <- rownames(tsneMat)
      }else{
        corDist <- as.dist(1-cor(mat4tsne_subset))
        tsneAUC <- Rtsne::Rtsne(corDist,
                                is_distance=TRUE,
                                perplexity=perplexity, ...)
        rownames(tsneAUC$Y) <- labels(corDist)
      }
      colnames(tsneAUC$Y) <- c("tsne1", "tsne2")
      tsneAUC$type <- paste0(aucType," ", nrow(mat4tsne_subset), " regulons (", dimsAsText, ", ", perplexity," perplexity)")
      
      fileName <- tsneFileName(filePrefix=filePrefix, aucType=aucType, nPcs=initial_dims, perpl=perplexity)
      if(aucType=="Binary")
      {
        tsneBinaryAUC <- tsneAUC
        saveRDS(tsneBinaryAUC, file=fileName)
      }else{
        saveRDS(tsneAUC, file=fileName)
      }
      return(fileName)
    }, error = function(e) {
      msg <- paste0("Couldn't create tSNE for: ", aucType,", ", dimsAsText,", ", perplexity, "perpl.\nError msg:", e)
      warning(msg)
      return(NULL)
    })
  }
  return(unlist(fileNames))
}

#' @export
tsneFileName <- function(scenicOptions=NULL, filePrefix=NULL, aucType=NULL, nPcs=NULL, perpl=NULL){
  if(is.null(filePrefix)) filePrefix <- getSettings(scenicOptions, "tSNE_filePrefix")
  if(is.null(aucType)) aucType <- getSettings(scenicOptions,"defaultTsne/aucType")
  if(is.null(nPcs)) nPcs <- getSettings(scenicOptions,"defaultTsne/dims")
  if(is.null(perpl)) perpl <- getSettings(scenicOptions,"defaultTsne/perpl")
  
  dimsAsText <- nPcs 
  if(tolower(nPcs)!="dist") dimsAsText <- paste0(sprintf("%02d", dimsAsText), "PCs")
  
  paste0(filePrefix, "_", aucType, "_", tolower(dimsAsText), "_", sprintf("%02d", perpl), "perpl.Rds")
}
