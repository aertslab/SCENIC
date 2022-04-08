#' @title runSCENIC_4_aucell_binarize
#' @description Step 4: Binarize the AUC (and, optional: re-cluster)
#' @param scenicOptions Fields used: TODO
#' @param skipBoxplot Whether to plot the boxplots
#' @param skipHeatmaps Whether to plot the Binary heatmaps
#' @param skipTsne Whether to calculate the binary t-SNE (and optionally, generate a HTML preview)
#' @param exprMat If !skipTsne, the expression matrix can be provided to plot the TF expression on the t-SNE (plotTsne_AUCellHtml)
#' @return The output is written in the folders 'int' and 'ouput'
#' @details See the detailed vignette explaining the internal steps.
#' @examples 
#' runSCENIC_4_aucell_binarize(scenicOptions)
#' @export
runSCENIC_4_aucell_binarize <- function(scenicOptions, 
                                        skipBoxplot=FALSE, skipHeatmaps=FALSE, skipTsne=FALSE, exprMat=NULL)
{
  nCores <- getSettings(scenicOptions, "nCores")
  regulonAUC <- tryCatch(loadInt(scenicOptions, "aucell_regulonAUC"),
                         error = function(e) {
                           if(getStatus(scenicOptions, asID=TRUE) < 3) 
                           e$message <- paste0("It seems the regulons have not been scored on the cells yet. Please, run runSCENIC_3_scoreCells() first.\n", 
                                               e$message)
                           stop(e)
  })
  thresholds <- loadInt(scenicOptions, "aucell_thresholds")
  thresholds <- getThresholdSelected(thresholds)
  
  # Assign cells
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(regulonAUC)[x,]>trh))
                                   }),names(thresholds))
  ### Convert to matrix (regulons with zero assigned cells are lost)
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"
  saveRDS(binaryRegulonActivity, file=getIntName(scenicOptions, "aucell_binary_full"))
  if(nrow(binaryRegulonActivity)==0) stop("No cells passed the binarization.")
  
  # Keep only non-duplicated thresholds
  # (e.g. only "extended" regulons if there is not a regulon based on direct annotation)
  binaryRegulonActivity_nonDupl <- binaryRegulonActivity[which(rownames(binaryRegulonActivity) %in% onlyNonDuplicatedExtended(rownames(binaryRegulonActivity))),]
  saveRDS(binaryRegulonActivity_nonDupl, file=getIntName(scenicOptions, "aucell_binary_nonDupl"))
  
  minCells <- ncol(binaryRegulonActivity) * .01
  msg <- paste0("Binary regulon activity: ",
                nrow(binaryRegulonActivity_nonDupl), " TF regulons x ",
                ncol(binaryRegulonActivity), " cells.\n(",
                nrow(binaryRegulonActivity), " regulons including 'extended' versions)\n",
                sum(rowSums(binaryRegulonActivity_nonDupl)>minCells),
                " regulons are active in more than 1% (", minCells, ") cells.")
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  
  if(!skipBoxplot)
  {
    .openDev(fileName=getOutName(scenicOptions, "s4_boxplotBinaryActivity"),
             devType=getSettings(scenicOptions, "devType"))
    par(mfrow=c(1,2))
    boxplot(rowSums(binaryRegulonActivity_nonDupl), main="nCells per regulon",
            sub='number of cells \nthat have the regulon active',
            col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)
    boxplot(colSums(binaryRegulonActivity_nonDupl), main="nRegulons per Cell",
            sub='number of regulons \nactive per cell',
            col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)
    dev.off()
  }
  
  ################################################################################
  # Binary activity heatmap
  if(!skipHeatmaps)
  {
    regulonSelection <- loadInt(scenicOptions, "aucell_regulonSelection", ifNotExists="null", verbose=FALSE)
    if(is.null(regulonSelection)) 
      regulonSelection <- regulonSelections(binaryRegulonActivity, binaryRegulonActivity_nonDupl, minCells)
    
    cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
    cellInfo <- data.frame(cellInfo)
    colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
    
    
    ### Plot heatmap:
    for(selRegs in names(regulonSelection$labels))
    {
      if(length(regulonSelection[[selRegs]])>0)
      {
        regulonSelection[[selRegs]] <- regulonSelection[[selRegs]][which(regulonSelection[[selRegs]] %in% rownames(binaryRegulonActivity))]
        binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]],,drop=FALSE]
        
        if(nrow(binaryMat)>0) 
        {
          fileName <- paste0(getOutName(scenicOptions, "s4_binaryActivityHeatmap"),selRegs)
          fileName <- .openDevHeatmap(fileName=fileName, devType=getSettings(scenicOptions, "devType"))
          
          rowv <- ifelse(nrow(binaryMat) >= 2, T, NA)
          colv <- ifelse(ncol(binaryMat) >= 2, T, NA)
          
          NMF::aheatmap(binaryMat, scale="none", revC=TRUE, main=selRegs,   
                        annCol=cellInfo[colnames(binaryMat),, drop=FALSE],
                        annColor=colVars,
                        Rowv=rowv,
                        Colv=colv,
                        color = c("white", "black"),
                        filename=fileName)
          if(getSettings(scenicOptions, "devType")!="pdf") dev.off()
        }else{
          if(getSettings(scenicOptions, "verbose")) message(paste0("No regulons to plot for regulon selection '", selRegs, "'. Skipping."))
        }
      }
    }
  }
  
  ################################################################################
  # Tsne - on binary activity
  if(!skipTsne)
  {
    tSNE_fileName <- tsneAUC(scenicOptions, aucType="Binary", filePrefix=getIntName(scenicOptions, "tsne_prefix"), onlyHighConf=FALSE) # default: nPcs, perpl, seed
    if(!is.null(tSNE_fileName))
    {
      tSNE <- readRDS(tSNE_fileName)
      
      # AUCell (activity) as html: 
      fileName <- getOutName(scenicOptions, "s4_binarytSNE_colAct")
      plotTsne_AUCellHtml(scenicOptions, exprMat, fileName, tSNE) #open the resulting html locally
      
      # Plot cell properties:
      sub <- ""; if("type" %in% names(tSNE)) sub <- paste0("t-SNE on ", tSNE$type)
      cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
      colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
      pdf(paste0(getOutName(scenicOptions, "s4_binarytSNE_colProps"),".pdf"))
      plotTsne_cellProps(tSNE$Y, cellInfo=cellInfo, colVars=colVars, cex=1, sub=sub)
      dev.off()
    }
  }
  
  # Finished. Update status.
  scenicOptions@status$current <- 4
  invisible(scenicOptions)
}


################################################################################
# Regulon orders/selection for plots
#' @export
regulonSelections <- function(binaryRegulonActivity, binaryRegulonActivity_nonDupl, minCells)
{
  #binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_full")
  #binaryRegulonActivity_nonDupl <- loadInt(scenicOptions, "aucell_binary_nonDupl")
  
  ### Select regulons:
  regulonSelection <- list(labels=c(all="All regulons \n (including duplicated regulons)",
                                    corr="Regulons with any other regulon correlated\n with abs(cor)>0.30 \n(and active in at least 1% of cells)",
                                    onePercent="Regulons active in more than 1% of cells",
                                    notCorr="Regulons with no other regulons correlated\n abs(cor)>0.30 \n or active in fewer than 1% of cells"))
  
  # All regulons.
  regulonSelection[["all"]] <- rownames(binaryRegulonActivity)
  
  # Active in > 1% cells
  regMinCells <- names(which(rowSums(binaryRegulonActivity_nonDupl) > minCells))
  regulonSelection[["onePercent"]] <- regMinCells
  
  # Correlation across regulons (based on binary cell activity)
  reguCor <- cor(t(binaryRegulonActivity_nonDupl[regMinCells,]))
  reguCor[which(is.na(reguCor))] <- 0
  diag(reguCor) <- 0
  
  # Regulons that co-ocurr in similar cells. If a regulon is relevant by itself it will not be shown, also check the regulons ignored.
  corrRegs <- names(which(rowSums(abs(reguCor) > 0.30) > 0))
  regulonSelection[["corr"]]  <- corrRegs
  missingRegs <- rownames(binaryRegulonActivity_nonDupl)[which(!rownames(binaryRegulonActivity_nonDupl) %in% corrRegs)]
  regulonSelection[["notCorr"]]  <- missingRegs
  saveRDS(regulonSelection, file=getIntName(scenicOptions, "aucell_regulonSelection"))
  
  ## Set regulon order (only plotting most correlated regulons)
  reguCor_dist <- as.dist(1-reguCor[corrRegs,corrRegs])
  if(length(reguCor_dist) >= 2) 
  {
    binaryRegulonOrder <- hclust(reguCor_dist)
    binaryRegulonOrder <- binaryRegulonOrder$labels[binaryRegulonOrder$order]
  } else 
  {
    binaryRegulonOrder <- labels(reguCor_dist)
  }
  saveRDS(binaryRegulonOrder, file=getIntName(scenicOptions, "aucell_binaryRegulonOrder"))
  
  return(regulonSelection)
}

