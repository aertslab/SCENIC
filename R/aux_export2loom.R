#' @title export2loom
#' @description Create .loom file with the results of the analysis to visualize them in http://scope.aertslab.org
#' @param scenicOptions scenicOptions object
#' @param dgem Digital expression matrix
#' @param hierarchy Labels for the hierarchy levels (up to three)
#' @return The .loom file (\code{file name indicated in "loomFile" slot in scenicOptions}) containing the following information:
#' \itemize{
#' \item dgem (digital expression matrix)
#' \item title: \code{getDatasetInfo(scenicOptions,"datasetTitle")}
#' \item genome:\code{ getDatasetInfo(scenicOptions,"org")}
#' \item regulons.AUC: \code{getAUC(loadInt(scenicOptions, "aucell_regulonAUC"))}
#' \item regulons: \code{loadInt(scenicOptions, "aucell_regulons")}
#' \item default.embedding: default t-SNE
#' \item default.embedding.name: default t-SNE type (e.g. AUC/binary and description)
#' \item column attributes: \code{getDatasetInfo(scenicOptions, "cellInfo")} and nGene (\code{colSums(dgem>0)}).
#' }
#' @seealso To add extra data (e.g. embeddings or clusters), see: \code{help(package="SCopeLoomR")}
#' @examples 
#' # devtools::install_github("aertslab/SCopeLoomR")
#' # scenicOptions@fileNames$output["loomFile",] <- "myAnalysis.loom"
#' # export2scope(scenicOptions, exprMat, hierarchy=c("SCENIC", "MouseBrain"))
#' @export 
export2loom <- export2scope <- function(scenicOptions, dgem, hierarchy=c("SCENIC", "", ""), 
                                        addAllTsnes=TRUE)
{
  regulonType="Motif" # TODO
  # TODO: what about if there is no dgem, but only normalized?
  
  if(length(hierarchy) > 3) stop("Maximum three hierarchy levels")
  if(length(hierarchy) < 3) hierarchy <- c(hierarchy, rep("", 5-length(hierarchy)))
  fileName <- getOutName(scenicOptions, "loomFile")
  if(file.exists(fileName)) stop("File '", fileName, "' already exists.")
  
  # TODO: ask max about order of samples tsne-expr-info
  suppressPackageStartupMessages(require(SCopeLoomR))

  # Default embedding (e.g. t-SNE or PCA coordinates)
  if(!file.exists(tsneFileName(scenicOptions))) stop(paste("Default 2D projection is not available:", tsneFileName(scenicOptions)))
  defaultTsne <- readRDS(tsneFileName(scenicOptions))
  defaultTsne_name <- paste("SCENIC t-SNE:", defaultTsne$type)
  defaultTsne <- defaultTsne$Y
  cellOrder <- rownames(defaultTsne)
  if(!all(colnames(dgem) == cellOrder)) warning("tSNE and matrix cell names do not match.")
  
  # Cell info:
  cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
  if(!all(rownames(cellInfo) == cellOrder)) warning("CellInfo cell names (order?) do not match.")
  for(cn in colnames(cellInfo)) if(is.numeric(cellInfo[,cn]) & is.character(cellInfo[,cn])) stop(paste0(cn, "should be either numeric or character/factor."))
  if(getSettings(scenicOptions, "verbose")) {
    message("The folowing cell metadata will be added:")
    print(data.frame(type=cbind(sapply(cellInfo, class))))
  }
  
  # Start adding...
  fileName <- getOutName(scenicOptions, "loomFile")
  loom <- build_loom(file.name=fileName,
                     dgem=dgem,
                     title=getDatasetInfo(scenicOptions,"datasetTitle"),
                     genome=getDatasetInfo(scenicOptions,"org"), # Just for user information, not used internally
                     default.embedding=defaultTsne,
                     default.embedding.name=defaultTsne_name)
  
  add_hierarchy(loom=loom, hierarchy=create_hierarchy(level.1.name=hierarchy[1], 
                                                      level.2.name=hierarchy[2], 
                                                      level.3.name=hierarchy[3]))
  
  # Known cell information/annotation  
  if(!is.null(cellInfo))
  {
    cellInfo <- data.frame(cellInfo)
    cellInfo <- cellInfo[,colnames(cellInfo) != "nGene", drop=FALSE]
    cellInfo <- cellInfo[,colnames(cellInfo) != "nUMI", drop=FALSE]
  
    # Add annotation
    for(cn in colnames(cellInfo))
    {
      isMetric <- FALSE
      isAnnotation <- FALSE
      if(is.character(cellInfo[,cn]) || is.factor(cellInfo[,cn])) 
      {
        isAnnotation <- TRUE
        cellInfo[,cn] <- as.character(cellInfo[,cn])
      }else{
        if(all(!is.na(as.numeric(as.character(cellInfo[,cn]))))) isMetric <- TRUE
      }
      add_col_attr(loom=loom, key=cn, value=cellInfo[,cn], as.annotation=isAnnotation, as.metric=isMetric)
    }
  }
  
  # Regulons AUC matrix
  regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
  if(!all(colnames(regulonAUC) == cellOrder)) warning("regulonAUC cell names (order?) do not match.")
  add_scenic_regulons_auc_matrix(loom=loom, regulons.AUC=AUCell::getAUC(regulonAUC), column.attr.name=paste0(regulonType, "RegulonsAUC"))
  
  # Regulons (gene list)
  regulons <- loadInt(scenicOptions, "aucell_regulons")
  motifEnrichment <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
  motifEnrichment <- motifEnrichment[grep("transfac_pro__", motifEnrichment[["motif"]], invert = T),] #TODO
  regulonThresholds <- loadInt(scenicOptions, "aucell_thresholds")
  if(!"aucThr" %in% names(regulonThresholds[[1]]))
  {
    # regulonThresholds <- lapply(regulonThresholds, function(x) list(selected=x))
    warning("The binarized regulon activity will not been added to the loom.")
    regulonThresholds <- NULL
  }
  add_scenic_regulons(loom=loom
                      , regulons=regulons
                      , regulon.threshold.assignments=regulonThresholds # Optional
                      , regulon.enrichment.table=motifEnrichment # Optional
                      , column.attr.name=paste0(regulonType, "Regulons")
                      )
  
  # # Alternative t-SNE #TODO
  # load(paste0(scenicDir, "int/3.3_tsneRegulonAUC_PCA_50pc50perpl.RData"))
  # add_embedding(loom=loom, embedding=tSNE$Y, name="SCENIC (t-SNE on AUC)")
  if(addAllTsnes)
  {
    tsnePath <- dirname(getSettings(scenicOptions, "tSNE_filePrefix"))
    allTsneFiles <- grep(pattern = paste0(basename(getSettings(scenicOptions, "tSNE_filePrefix")), ".*\\.Rds")
         , x=list.files(tsnePath)
         , value = T)
    
    for(tsneFile in allTsneFiles)
    {
      tryCatch(
      add_embedding(loom=loom, embedding=readRDS(file.path(tsnePath, tsneFile))$Y, name=paste0("SCENIC: ",  gsub(".Rds","", tsneFile, fixed=TRUE)))
      , error=function(e) print(paste0("Cannot add tsne: ", e$message)))
    }
  }
  
  finalize(loom=loom)
  if(getSettings(scenicOptions, "verbose")) message("Loom file saved as:\t", fileName)
}
