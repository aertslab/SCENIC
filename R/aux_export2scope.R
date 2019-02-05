#' @title export2scope
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
export2scope <- function(scenicOptions, dgem, hierarchy=c("SCENIC", "", ""), addAllTsnes=TRUE)
{
  # TODO: what about if there is no dgem, but only normalized?
  
  if(length(hierarchy) > 3) stop("Maximum three hierarchy levels")
  if(length(hierarchy) < 3) hierarchy <- c(hierarchy, rep("", 5-length(hierarchy)))
  
  # TODO: ask max about order of samples tsne-expr-info
  suppressPackageStartupMessages(require(SCopeLoomR))

  # Default embedding (e.g. t-SNE or PCA coordinates)
  defaultTsne <- readRDS(tsneFileName(scenicOptions))
  defaultTsne_name <- paste("SCENIC t-SNE:", defaultTsne$type)
  defaultTsne <- defaultTsne$Y
  
  fileName <- getOutName(scenicOptions, "loomFile")
  if(file.exists(fileName)) stop("File '", getOutName(scenicOptions, "loomFile"), "' already exists.")
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
  cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
  if(!is.null(cellInfo))
  {
    cellInfo <- data.frame(cellInfo)
    cellInfo <- cellInfo[,colnames(cellInfo) != "nGene", drop=FALSE]
    cellInfo <- cellInfo[,colnames(cellInfo) != "nUMI", drop=FALSE]
  
    # Add annotation
    for(cn in colnames(cellInfo))
    {
      add_col_attr(loom=loom, key=cn, value=cellInfo[,cn])
    }
  }
  
  # Regulons AUC matrix
  regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
  add_scenic_regulons_auc_matrix(loom=loom, regulons.AUC=AUCell::getAUC(regulonAUC))
  
  # Regulons (gene list)
  regulons <- loadInt(scenicOptions, "aucell_regulons")
  motifEnrichment <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
  motifEnrichment <- motifEnrichment[grep("transfac_pro__", motifEnrichment[["motif"]], invert = T),] #TODO
  regulonThresholds <- loadInt(scenicOptions, "aucell_thresholds")
  add_scenic_regulons(loom=loom
                      , regulons=regulons
                      , regulon.threshold.assignments=regulonThresholds # Optional
                      , regulon.enrichment.table=motifEnrichment # Optional
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
      add_embedding(loom=loom, embedding=readRDS(file.path(tsnePath, tsneFile))$Y, name=paste0("SCENIC: ", file.path(tsneFile)))
      , error=function(e) print(paste0("Cannot add tsne: ", e$message)))
    }
  }
  
  finalize(loom=loom)
}
