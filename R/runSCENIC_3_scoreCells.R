################################################################################
# Step 3. Analyzing the network activity in each individual cell
################################################################################
#' @title runSCENIC_3_scoreCells
#' @description Step 3: AUCell (scoring the regulons on the individual cells) 
#' @param scenicOptions Fields used: TODO
#' @param exprMat Expression matrix
#' @param skipBinaryThresholds Whether to skip the automatic binarization step
#' @param skipHeatmap Whether to plot the AUC heatmap
#' @param skipTsne Whether to plot the t-SNE
#' @return The output is written in the folders 'int' and 'ouput'
#' @details See the detailed vignette explaining the internal steps.
#' @examples 
#' runSCENIC_3_scoreCells(scenicOptions)
#' @export
runSCENIC_3_scoreCells <- function(scenicOptions, exprMat, 
                               skipBinaryThresholds=FALSE, skipHeatmap=FALSE, skipTsne=FALSE)
{
  nCores <- getSettings(scenicOptions, "nCores")
  
  ################################################################
  ## Prepare regulons
  regulons <- loadInt(scenicOptions, "regulons")
  regulons <- regulons[order(lengths(regulons), decreasing=TRUE)]
  regulons <- regulons[lengths(regulons)>=10]
  if(length(regulons) <2)  stop("Not enough regulons with at least 10 genes.")
  
  # Add the TF to the regulon (keeping it only once) & rename regulon
  regulons <- setNames(lapply(names(regulons), function(tf) sort(unique(c(gsub("_extended", "", tf), regulons[[tf]])))), names(regulons))
  names(regulons) <- paste(names(regulons), " (",lengths(regulons), "g)", sep="")
  saveRDS(regulons, file=getIntName(scenicOptions, "aucell_regulons"))
  
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tStep 3. Analyzing the network activity in each individual cell")
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  msg <- paste0("\nNumber of regulons to evaluate on cells: ", length(regulons),
                "\nBiggest (non-extended) regulons: \n",
                paste("\t", grep("_extended",names(regulons),invert = T, value = T)[1:10], collapse="\n")) # TODO maxlen?
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  ################################################################
  # AUCell
  library(AUCell)
  # 1. Create rankings
  set.seed(getSettings(scenicOptions,"seed"))
  .openDev(fileName=getIntName(scenicOptions, "aucell_genesStatsPlot"), 
          devType=getSettings(scenicOptions, "devType"))
    aucellRankings <- AUCell_buildRankings(exprMat, nCores=nCores, 
                          plotStats=TRUE, verbose=getSettings(scenicOptions, "verbose"))
    abline(v=aucellRankings@nGenesDetected["1%"], col="skyblue3", lwd=5, lty=3) #TODO check if similar to previous version
  dev.off()
  saveRDS(aucellRankings, file=getIntName(scenicOptions, "aucell_rankings"))
    
  # 2. Calculate AUC
  regulonAUC <- AUCell_calcAUC(regulons, aucellRankings, 
              aucMaxRank=aucellRankings@nGenesDetected["1%"], nCores=nCores)
    
  # Order the modules by similarity, for easier exploration in the upcoming steps & save
  variableRegulons <- names(which(apply(getAUC(regulonAUC), 1, sd) > 0))
  reguDist <-as.dist(1-cor(t(getAUC(regulonAUC)[variableRegulons,]), method="spear"))
  reguClust <- hclust(reguDist, method="ward.D2")
  regulonClusters <- setNames(dynamicTreeCut::cutreeDynamic(reguClust, distM=as.matrix(reguDist), verbose = FALSE), reguClust$labels)
  regulonOrder <- reguClust$labels[reguClust$order]
  regulonOrder <- regulonOrder[order(regulonClusters[regulonOrder], decreasing = TRUE)]
  
  regulonAUC <- regulonAUC[regulonOrder,]
  saveRDS(regulonAUC, file=getIntName(scenicOptions, "aucell_regulonAUC"))
  
  # 3. Default thresholds
  cells_AUCellThresholds <- NULL
  if(!skipBinaryThresholds)
  {
    cells_AUCellThresholds <- AUCell_exploreThresholds(regulonAUC, 
                          smallestPopPercent=getSettings(scenicOptions,"aucell/smallestPopPercent"),
                          assignCells=TRUE, plotHist=FALSE, 
                          verbose=FALSE, nCores=nCores)
    saveRDS(cells_AUCellThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
    
    # Get cells assigned to each regulon
    regulonsCells <- getAssignments(cells_AUCellThresholds)
    
    ### Save threshold info as text (e.g. to edit/modify...)
    trhAssignment <- getThresholdSelected(cells_AUCellThresholds)
    trhAssignment <- signif(trhAssignment,3)
    commentsThresholds <- sapply(cells_AUCellThresholds, function(x) unname(x$aucThr$comment))
    
    table2edit <- cbind(regulon=names(trhAssignment),
                        threshold=trhAssignment,
                        nCellsAssigned=lengths(regulonsCells)[names(trhAssignment)],
                        AUCellComment=commentsThresholds,
                        nGenes=gsub("[\\(g\\)]", "", regmatches(names(cells_AUCellThresholds), gregexpr("\\(.*?\\)", names(cells_AUCellThresholds)))),
                        clusteringOrder=1:length(trhAssignment),
                        clusterGroup=regulonClusters[names(trhAssignment)],
                        onlyNonDuplicatedExtended=(names(trhAssignment) %in% onlyNonDuplicatedExtended(names(trhAssignment))),
                        personalNotes="")
    write.table(table2edit, file=getIntName(scenicOptions, "aucell_thresholdsTxt"), row.names=F, quote=F, sep="\t")
    rm(trhAssignment)
  }
  ####################################################################
  # Plots
  if(!skipHeatmap){
    nCellsHeatmap <- min(500, ncol(regulonAUC))
    cells2plot <- sample(colnames(regulonAUC), nCellsHeatmap)
    
    cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")   #TODO check if exists, if not... create/ignore?
    if(!is.null(cellInfo)) cellInfo <- data.frame(cellInfo)[cells2plot,,drop=F]
    colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
    
    fileName <- getOutName(scenicOptions, "s3_AUCheatmap")
    
    fileName <- .openDevHeatmap(fileName=fileName, devType=getSettings(scenicOptions, "devType"))
    NMF::aheatmap(getAUC(regulonAUC)[,cells2plot],
                  annCol=cellInfo,
                  annColor=colVars,
                  main="AUC",
                  sub=paste("Subset of",nCellsHeatmap," random cells"),
                  filename=fileName)
    .closeDevHeatmap(devType=getSettings(scenicOptions, "devType"))
  }
  ####################################################################
  # Plots
  if(!skipTsne){
    tSNE_fileName <- tsneAUC(scenicOptions, aucType="AUC", onlyHighConf=FALSE) # default: nPcs, perpl, seed, tsne prefix
    tSNE <- readRDS(tSNE_fileName)
    
    # AUCell (activity) plots with the default tsne, as html: 
    fileName <- getOutName(scenicOptions, "s3_AUCtSNE_colAct")
    plotTsne_regulonActivityHTML(scenicOptions, exprMat, fileName, tSNE) #open the resulting html locally

    # Plot cell properties:
    sub <- ""; if("type" %in% names(tSNE)) sub <- paste0("t-SNE on ", tSNE$type)
    cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null") 
    colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
    pdf(paste0(getOutName(scenicOptions, "s3_AUCtSNE_colProps"),".pdf"))
    plotTsne_cellProps(tSNE$Y, cellInfo=cellInfo, colVars=colVars, cex=1, sub=sub)
    dev.off()
  }
}

