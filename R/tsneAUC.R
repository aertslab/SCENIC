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
#' #todo
#' # tSNE <- readRDS(defaultTsneFileName(scenicOptions))
#' @export
tsneAUC <- function(scenicOptions, aucType=NULL, nPcs=NULL, perpl=NULL, filePrefix="tSNE", seed=NULL, onlyHighConf=FALSE, ...)
{
  if(is.null(aucType)) aucType <- getSettings(scenicOptions,"defaultTsne/aucType")
  if(is.null(nPcs)) nPcs <- getSettings(scenicOptions,"defaultTsne/dims")
  if(is.null(perpl)) perpl <- getSettings(scenicOptions,"defaultTsne/dims")
  
  if(is.null(seed)) seed <- getSettings(scenicOptions,"seed")
  
  if(!aucType %in% c("Binary", "AUC")) stop('aucType should be either "Binary" or "AUC"')
    
  if(aucType=="Binary") mat4tsne <- loadInt(scenicOptions, "aucell_binary_nonDupl") 
  if(aucType=="AUC") mat4tsne <- getAUC(loadInt(scenicOptions, "aucell_regulonAUC"))

  msg <- paste0(format(Sys.time(), "%H:%M"), "\tCalculating t-SNE")
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  if(onlyHighConf)
  {
    mat4tsne_subset <- mat4tsne[grep("_extended",rownames(mat4tsne),invert=T, value=T),] 
  }else{
    mat4tsne_subset <- mat4tsne[onlyNonDuplicatedExtended(rownames(mat4tsne)),]
  }
  
  fileNames <- c()
  for(perplexity in perpl)
  {
    for(initial_dims in nPcs)
    {
      # PCA-based t-SNE
      set.seed(seed)
      if(aucType=="Binary") mat4tsne_subset <- jitter(mat4tsne_subset, factor=1)
      tsneMat <- unique(t(mat4tsne_subset))  
      tsneAUC <- Rtsne::Rtsne(tsneMat,
                              initial_dims=initial_dims,
                              perplexity=perplexity, ...) #TODO tSNE function?
      tsneAUC$type <- paste0(aucType," ", nrow(mat4tsne_subset), " regulons (", initial_dims, " PCs, ", perplexity," perplexity)")
      rownames(tsneAUC$Y) <- rownames(tsneMat)
      colnames(tsneAUC$Y) <- c("tsne1", "tsne2")
      
      fileName <- paste0(filePrefix, "_",aucType, "_", initial_dims, "pcs_", perplexity, "perpl.Rds") #TODO
      if(aucType=="Binary") 
      {
        tsneBinaryAUC <- tsneAUC
        saveRDS(tsneBinaryAUC, file=fileName)
      }else{
        saveRDS(tsneAUC, file=fileName)
      }
      fileNames <- c(fileNames, fileName)
    }
  }
  
  return(fileNames)
}

#' @export
defaultTsneFileName <- function(scenicOptions){
  filePrefix <- getIntName(scenicOptions, "tsne_prefix")
  aucType <- getSettings(scenicOptions,"defaultTsne/aucType")
  nPcs <- getSettings(scenicOptions,"defaultTsne/dims")
  perpl <- getSettings(scenicOptions,"defaultTsne/dims")
  
  paste0(filePrefix, "_",aucType, "_", nPcs, "pcs_", perpl, "perpl.Rds") #TODO
}
