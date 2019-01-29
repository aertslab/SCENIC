
################################################################################
# Step 1. Part 2: Creating TF modules (potential TF-targets)
################################################################################


#' @title runSCENIC_1_coexNetwork2modules
#' @description Step 1: Convert the output from GENIE3/GRNBoost to co-expression modules
#' @param scenicOptions Fields used: TODO
#' @return The output is written in the folders 'int' and 'ouput'
#' @details See the detailed vignette explaining the internal steps.
#' @examples 
#' runSCENIC_1_coexNetwork2modules(scenicOptions)
#' @export
runSCENIC_1_coexNetwork2modules <- function(scenicOptions)
{
  linkList <- loadInt(scenicOptions, "genie3ll")
  if(!all(colnames(linkList) == c("TF", "Target", "weight"))) 
    stop('The link list colnames should be "TF", "Target", "weight"')
  
  uniquePairs <- nrow(unique(linkList[,c("TF", "Target")]))
  if(uniquePairs < nrow(linkList)) 
    stop("There are duplicated regulator-target (gene id/name) pairs in the input link list.")
  
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tCreating TF modules")
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  print(quantile(linkList$weight, probs=c(0.75, 0.90)))
  .openDev(fileName=getIntName(scenicOptions, "genie3weighPlot"), 
           devType=getSettings(scenicOptions, "devType"))
    plot(linkList$weight[1:1000000], type="l", ylim=c(0, max(linkList$weight)), main="Weight of the links",
         ylab="Weight", xlab="Links sorted decreasingly")
    abline(h=0.001, col="blue") # Threshold
    #sum(linkList$weight>0.001)/nrow(linkList)
  dev.off()
  
  # Keep only genes with weight > threshold
  linkList_001 <- linkList[which(linkList[,"weight"]>getSettings(scenicOptions, "modules/weightThreshold")),]
  if(getSettings(scenicOptions, "verbose")) message("Number of links between TFs and targets: ", nrow(linkList_001))
  
  #### Create the gene-sets & save:
  tfModules <- list()
  
  linkList_001$TF <- as.character(linkList_001$TF)
  linkList_001$Target <- as.character(linkList_001$Target)
  
  ### Create TF-modules:
  # 1: Weight > 0.001 (filtered in previous step)
  tfModules[["w001"]] <- split(linkList_001$Target, factor(linkList_001$TF))
  
  # 2: Weight > 0.005
  llminW <- linkList_001[which(linkList_001[,"weight"]>0.005),]
  tfModules[["w005"]] <- split(llminW$Target, factor(llminW$TF))
  
  # 3: Top 50 targets for each TF
  # ("w001" should be ordered decreasingly by weight)
  tfModules[["top50"]] <- lapply(tfModules[["w001"]], function(x) x[1:(min(length(x), 50))])
  
  # 4-6: Top regulators per target
  # (linkList_001 should be ordered by weight!)
  linkList_001_byTarget <- split(linkList_001, factor(linkList_001$Target))
  
  nTopTfs <- c(5, 10, 50)
  nTopTfs <- setNames(nTopTfs, paste("top", nTopTfs, "perTarget", sep=""))
  
  #library(reshape2); library(data.table)
  topTFsperTarget <- lapply(linkList_001_byTarget, function(llt) {
    nTFs <- nTopTfs[which(nTopTfs <= nrow(llt))]
    reshape2::melt(lapply(nTFs, function(x) llt[1:x,"TF"]))
  })
  
  topTFsperTarget <- topTFsperTarget[which(!sapply(sapply(topTFsperTarget, nrow), is.null))]
  topTFsperTarget.asDf <-  data.frame(data.table::rbindlist(topTFsperTarget, idcol=TRUE))
  # topTFsperTarget.asDf <- apply(topTFsperTarget.asDf, 2, as.character)
  colnames(topTFsperTarget.asDf) <- c("Target", "TF", "method")
  
  # Merge the all the gene-sets:
  tfModules.melted <- reshape2::melt(tfModules)
  colnames(tfModules.melted) <- c("Target", "TF", "method")
  tfModules <- rbind(tfModules.melted, topTFsperTarget.asDf)
  rm(tfModules.melted); rm(topTFsperTarget.asDf)
  tfModules$TF <- as.character(tfModules$TF)
  tfModules$Target <- as.character(tfModules$Target)
  
  # Basic counts:  #TODO add comment
  if(getSettings(scenicOptions, "verbose")) 
    print(
      rbind(nTFs=length(unique(tfModules$TF)),
            nTargets=length(unique(tfModules$Target)),
            nGeneSets=nrow(unique(tfModules[,c("TF","method")])),
            nLinks=nrow(tfModules))
    )
  
  ### Add correlation to split into positive- and negative-correlated targets
  corrMat <- loadInt(scenicOptions, "corrMat")
  # Keep only correlation between TFs and potential targets
  tfs <- unique(tfModules$TF)
  missingTFs <- tfs[which(!tfs %in% rownames(corrMat))]
  if(length(missingTFs) >0 ) 
  { 
    warning("The following TFs are missing from the correlation matrix: ", paste(missingTFs, collapse=", "))
    
    tfs <- tfs[which(tfs %in% rownames(corrMat))]
    corrMat <- corrMat[tfs,]
  }
  
  # Add correlation to the table
  # "corr" column: 1 if the correlation between the TF and the target is > 0.03, -1 if the correlation is < -0.03 and 0 otherwise.
  tfModules_byTF <- split(tfModules, as.factor(tfModules$TF))
  tfModules_withCorr_byTF <- lapply(tfModules_byTF[tfs[1:4]], function(tfGeneSets)
  {
    tf <- as.character(unique(tfGeneSets$TF))
    targets <- as.character(tfGeneSets$Target)
    cbind(tfGeneSets, corr=c(as.numeric(corrMat[tf,targets] > 0.03) - as.numeric(corrMat[tf,targets] < -0.03)))
  })
  tfModules_withCorr <- data.frame(data.table::rbindlist(tfModules_withCorr_byTF))
  if(length(missingTFs) >0 )
  { 
    tfModules_withCorr <- rbind(tfModules_withCorr, data.frame(tfModules[tfModules$TF %in% missingTFs,], corr=NA))
  }
  saveRDS(tfModules_withCorr, file=getIntName(scenicOptions, "tfModules_asDF"))
}
