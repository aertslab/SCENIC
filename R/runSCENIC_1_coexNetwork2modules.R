
################################################################################
# Step 1. Part 2: Creating TF modules (potential TF-targets)
################################################################################


#' @title runSCENIC_1_coexNetwork2modules
#' @description Step 1: Convert the output from GENIE3/GRNBoost to co-expression modules
#' @param scenicOptions Fields used: TODO
#' @param weightThreshold Global weight threshold. Links with lower weight will be ignored in ALL modules.
#' @return The output is written in the folders 'int' and 'ouput'
#' @details See the detailed vignette explaining the internal steps.
#' @examples 
#' runSCENIC_1_coexNetwork2modules(scenicOptions)
#' @export
runSCENIC_1_coexNetwork2modules <- function(scenicOptions, 
                                            weightThreshold=getSettings(scenicOptions, "modules/weightThreshold"),
                                            topThr=c(0.005), # or: quantile([,weightCol], c(.75, .90))
                                            nTopTfs=c(5, 10, 50),
                                            nTopTargets=50, 
                                            aFun=topPerTf, 
                                            corrThr=0.03,
                                            linkList=NULL,
                                            corrMat=NULL,
                                            weightCol="weight",
                                            verbose=getSettings(scenicOptions, "verbose"))
{
  if(is.null(linkList))
  {
      linkList <- loadInt(scenicOptions, "genie3ll")
  }
  if(!all(c("TF", "Target", weightCol) %in% colnames(linkList)))
    stop('The link list colnames should be "TF", "Target", "', weightCol,'"')
  
  cntPairs <- table(table(linkList[,"TF"],linkList[,"Target"]))
  if(any(names(cntPairs)>1))
    stop("There are duplicated regulator-target (gene id/name) pairs in the input link list.")
  
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tCreating TF modules")
  if(verbose) message(msg)
  
  .openDev(fileName=getIntName(scenicOptions, "genie3weighPlot"), 
           devType=getSettings(scenicOptions, "devType"))
    plot(linkList[,weightCol][1:min(1000000,nrow(linkList))], type="l", ylim=c(0, max(linkList[,weightCol])), main="Weight of the links",
         ylab=weightCol, xlab="Links sorted decreasingly")
    abline(h=0.001, col="blue") # Threshold
    #sum(linkList[,weightCol]>0.001)/nrow(linkList)
  dev.off()
  
  # Keep only genes with weight > threshold
  linkList <- linkList[which(linkList[,weightCol] >= weightThreshold),]
  # weightStats <- quantile(linkList[,weightCol], probs=)
  if(verbose) 
  {
    # print(weightStats)
    message("Number of links between TFs and targets (weight>=",weightThreshold,"): ", nrow(linkList))
  }
  
  #### Create the gene-sets & save:
  tfModules <- list()
  linkList$TF <- as.character(linkList$TF)
  linkList$Target <- as.character(linkList$Target)
  
  ### Create TF-modules:
  # 1: Weight > 0.001 (filtered in previous step)
  allName <- paste0("w", format(weightThreshold, scientific=FALSE))
  tfModules[[allName]] <- split(linkList$Target, factor(linkList$TF))
  
  # 2: Weight > 0.005
  if(!is.null(topThr))
  {
    topThr <- setNames(topThr, paste("w", format(topThr, scientific=FALSE), sep=""))
    for(i in seq_along(topThr))
    {
      llminW <- linkList[which(linkList[,weightCol]>topThr[i]),]
      tfModules[[names(topThr)[i]]] <- split(llminW$Target, factor(llminW$TF))
    }
  }
    
  # 3: Top XX targets for each TF
  # (tfModules[[allName]] should be ordered decreasingly by weight)
  if(!is.null(nTopTargets))
  {
    nTopTargets <- setNames(nTopTargets, paste("top", nTopTargets, sep=""))
    for(i in seq_along(nTopTargets))
    {
      tfModules[[names(nTopTargets)[i]]] <- lapply(tfModules[[allName]], function(x) x[1:(min(length(x), nTopTargets[i]))]) 
    }
  }
  tfModules.melted <- reshape2::melt(tfModules)
  colnames(tfModules.melted) <- c("Target", "TF", "method")
  
  # 4-6: Top regulators per target
  # (linkList should be ordered by weight!)
  topTFsperTarget.asDf <- NULL
  if(!is.null(nTopTfs))
  {
    linkList_byTarget <- split(linkList, factor(linkList$Target))
    
    nTopTfs <- setNames(nTopTfs, paste("top", nTopTfs, "perTarget", sep=""))
    topTFsperTarget <- lapply(linkList_byTarget, function(llt) {
      nTFs <- nTopTfs[which(nTopTfs <= nrow(llt))]
      reshape2::melt(lapply(nTFs, function(x) llt[1:x,"TF"]))
    })
    
    topTFsperTarget <- topTFsperTarget[which(!sapply(sapply(topTFsperTarget, nrow), is.null))]
    topTFsperTarget.asDf <-  data.frame(data.table::rbindlist(topTFsperTarget, idcol=TRUE))
    # topTFsperTarget.asDf <- apply(topTFsperTarget.asDf, 2, as.character)
    colnames(topTFsperTarget.asDf) <- c("Target", "TF", "method")
  }
  
  # 7: Default function replaces "3" in most cases...
  byFun <- NULL
  if(!is.null(aFun))
  {
    byFun <- aFun(linkList, weightCol) 
  }
  
  # Merge the all the gene-sets:
  tfModules <- rbind(tfModules.melted, topTFsperTarget.asDf, byFun)
  rm(tfModules.melted); rm(topTFsperTarget.asDf)
  tfModules$TF <- as.character(tfModules$TF)
  tfModules$Target <- as.character(tfModules$Target)
  
  # Basic counts:  #TODO add comment
  if(verbose) 
    print(
      rbind(nTFs=length(unique(tfModules$TF)),
            nTargets=length(unique(tfModules$Target)),
            nGeneSets=nrow(unique(tfModules[,c("TF","method")])),
            nLinks=nrow(tfModules))
    )
  
  ### Add correlation to split into positive- and negative-correlated targets
  if(is.null(corrMat))
  {
      corrMat <- loadFile(scenicOptions, getIntName(scenicOptions, "corrMat"), verbose=FALSE, ifNotExists="null")
  }
  
  if(!is.null(corrMat))
  {
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
    tfModules_withCorr_byTF <- lapply(tfModules_byTF[tfs], function(tfGeneSets)
    {
      tf <- as.character(unique(tfGeneSets$TF))
      targets <- as.character(tfGeneSets$Target)
      cbind(tfGeneSets, corr=c(as.numeric(corrMat[tf,targets] > corrThr) - as.numeric(corrMat[tf,targets] < -corrThr)))
    })
    tfModules_withCorr_byTF <- tfModules_withCorr_byTF[which(lengths(tfModules_withCorr_byTF)>0)]
    tfModules_withCorr <- data.frame(data.table::rbindlist(tfModules_withCorr_byTF))
    if(length(missingTFs) >0 )
    { 
      tfModules_withCorr <- rbind(tfModules_withCorr, data.frame(tfModules[tfModules$TF %in% missingTFs,], corr=NA))
    }
  }else{
    tfModules_withCorr <-  data.frame(tfModules, corr=NA)
    if(verbose) message("Correlation information not available. It will not be added to the modules.")
  }
  
  saveRDS(tfModules_withCorr, file=getIntName(scenicOptions, "tfModules_asDF"))
  
  # Finished. Update status.
  scenicOptions@status$current <- 1
  invisible(scenicOptions)
}

#' @export
topPerTf <- function(ll, weightCol)
{
  ll <- split(ll, factor(as.character(ll$TF)))
  tptf <- setNames(lapply(names(ll), function(tf){
    tfMean <- mean(ll[[tf]][,weightCol])
    tfSd <- sd(ll[[tf]][,weightCol])
    list(top1sd=ll[[tf]][which(ll[[tf]][,weightCol] >= tfMean+tfSd),"Target"],
         top3sd=ll[[tf]][which(ll[[tf]][,weightCol] >= tfMean+(3*tfSd)),"Target"])
    # names(tmp) <- paste0(tf, "_", names(tmp))
    # tmp
  }),names(ll))
  tptf <- reshape2::melt(tptf)
  colnames(tptf)[which(colnames(tptf)=="value")] <- "Target"
  colnames(tptf)[which(colnames(tptf)=="L1")] <- "TF"
  colnames(tptf)[which(colnames(tptf)=="L2")] <- "method"
  tptf <- tptf[,c("Target", "TF", "method")]
  return(tptf)
}

# simplifyedMods <- mergeOverlappingModules(tfModules_withCorr)
# table(simplifyedMods$TF, simplifyedMods$method)

#' @export
mergeOverlappingModules <- function(tfModulesDf, minJakkardInd=.8)
{
  if(minJakkardInd<=0) stop("minJakkardInd should be > 0")
  tfModulesDf <- split(tfModulesDf, tfModulesDf$TF)
  
  # st <- list()
  for(tf in names(tfModulesDf))
  {
    tmp <- tfModulesDf[[tf]]
    tmp <- split(tmp$Target, tmp$method)
    st <- matrix(0, ncol=length(names(tmp)), nrow=length(names(tmp)), dimnames=list(names(tmp), names(tmp)))
    for(i in names(tmp))
      for(j in names(tmp))
      {
        if(i>j)
          st[i,j] <- signif(length(intersect(tmp[[i]], tmp[[j]]))/length(unique(c(tmp[[i]], tmp[[j]]))),2)
      }
    
    toMerge <- which(st>=minJakkardInd,arr.ind=T)
    if(nrow(toMerge)>0)
    {
      for(k in order(-st[toMerge])) 
      {
        mOne <- names(tmp)[toMerge[k,"row"]]
        mTwo <- names(tmp)[toMerge[k,"col"]] 
        newName <- paste0(mOne,"And",mTwo)
        if(all(c(mOne,mTwo) %in% tfModulesDf[[tf]]$method)) tfModulesDf[[tf]][which(tfModulesDf[[tf]]$method %in% c(mOne,mTwo)),"method"] <- newName
      }
      tfModulesDf[[tf]] <- unique(tfModulesDf[[tf]])
    } 
  }
  tfModulesDf <- data.frame(data.table::rbindlist(tfModulesDf))
  return(tfModulesDf)
}



