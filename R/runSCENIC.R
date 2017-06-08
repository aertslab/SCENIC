# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @import GENIE3
#' @import AUCell
#' @import RcisTarget
#' @import Biobase
#' @import BiocGenerics
#' @import data.table
#' @import reshape2
#'
#' @title Run SCENIC (steps 1.2 trhough 4)
#' @description Runs SCENIC steps after GENIE3.
#' This function should be equivalent to running the following tutorials sequentially :
#' - Step1.2_Regulons.Rmd:  Identifying regulons (direct TF targets) based on DNA motif enrichment.
#' - Step3.1_NwActivity.Rmd: Analyzing the network activity in each individual cell (part 1: calculate AUC).
#' - Step3.2_BinaryNwActivity.Rmd: Creating the binary activity matrix (convert the network activity into ON/OFF).
#' - Step4_Clustering.Rmd: Identify stable cell states based on their gene regulatory network activity (cell clustering)
#' (e.g. all steps and parameters are the same)
#' @param exprMat expression matrix
#' @param org "mm9" for mouse, or "hg19" for human. Determines the databases used for RcisTarget
#' @param cellInfo Phenodata information for the cells
#' @param colVars Colors to use to plot the cellInfo (formatted as a list, same format as in NMF::aheatmap)
#' @param stepsToRun Determines the Steps of SCENIC to run (matches the numbers/order in the title of the tutorials). By default: c("1.2", "2", "3.1", "3.2", "4").
#' @param nCores Number of cores to use for computation
#' @param seed Seed to use for "random" computations (for reproducibility of the results)
#' @param verbose TRUE/FALSE. Whether to show progress messages on screen. In either case, the messages will be saved in the log file: output/runSCENIC_log.txt
#' @param showAlternativeTsnes TRUE/FALSE. Whether to plot t-SNEs with diferent parameters.
#' In general it is recommended to plot multiple t-SNEs to evaluate the stability of the results, but it might take long time to run. For big datasets, it might be better to manually choose/try the parameters following the individual tutorials.
#' @return No value is returned. The outputs of this function are saved in the folder "output" the current working directory. Other files are saved in the folder "int" (e.g. for further exploration or to perform alternative analyses).
#' @details ...
#' @seealso List of vignettes included in the package: \code{vignette(package="SCENIC")}
#' @example inst/examples/runSCENIC.R
#' @export
runSCENIC <- function(exprMat=NULL, org=NULL, cellInfo=NULL, colVars=NULL,
                      stepsToRun=c("1.2", "2", "3.1", "3.2", "4"),
                      nCores=4, seed=123, verbose=TRUE, showAlternativeTsnes=TRUE)
{
  # Load libraries
  # suppressPackageStartupMessages({
  #   library(Biobase)
  #   library(BiocGenerics)
  #   library(data.table)
  #   library(reshape2)
  #
  #   library(GENIE3)
  #   library(AUCell)
  #   library(RcisTarget)
  # })

  # Check existing files
  dir.create("output", showWarnings=FALSE)
  sink("output/runSCENIC_log.txt", split=verbose)

  if(!file.exists("int")){
    stop("The output from GENIE3 should be located in a folder named 'int' in the current working directory.")
  }
  if("1.2" %in% stepsToRun)
  {
    if(!file.exists("int/1.3_GENIE3_weightMatrix.RData") && !file.exists("int/1.5_GENIE3_linkList.RData")){
      stop("The output from GENIE3 should be stored in 'int/1.3_GENIE3_weightMatrix.RData' or 'int/1.5_GENIE3_linkList.RData'")
    }
    if(!file.exists("int/1.4_corrMat.RData")){
      stop("The output from the correlation should be stored in 'int/1.4_corrMat.RData'")
    }
  }

  stepsToRun <- as.character(stepsToRun)
  ## Motif databases
  if(!org %in% c("mm9", "hg19", NULL)) stop("Organism not valid. Choose either 'hg19' (human), 'mm9' (mouse) or NULL.")
  if((org=="mm9") & ("2" %in% stepsToRun))
  {
    message(format(Sys.time(), "%H:%M"), "\tLoading mouse (mm9) databases.")
    library(RcisTarget.mm9.motifDatabases.20k)

    # Motif rankings (genes x motifs)
    data(mm9_500bpUpstream_motifRanking)
    data(mm9_10kbpAroundTss_motifRanking)
    motifRankings <- list()
    motifRankings[["500bp"]] <- mm9_500bpUpstream_motifRanking
    motifRankings[["10kbp"]] <- mm9_10kbpAroundTss_motifRanking

    # Motif annotation (TFs)
    data(mm9_direct_motifAnnotation)
    direct_motifAnnotation <- mm9_direct_motifAnnotation
    data(mm9_inferred_motifAnnotation) # optional
    inferred_motifAnnotation <- mm9_inferred_motifAnnotation
  }

  if((org=="hg19") & ("2" %in% stepsToRun))
  {
    message(format(Sys.time(), "%H:%M"), "\tLoading human (hg19) databases.")
    library(RcisTarget.hg19.motifDatabases.20k)

    # Motif rankings (genes x motifs)
    data(hg19_500bpUpstream_motifRanking)
    data(hg19_10kbpAroundTss_motifRanking)
    motifRankings <- list()
    motifRankings[["500bp"]] <- hg19_500bpUpstream_motifRanking
    motifRankings[["10kbp"]] <- hg19_10kbpAroundTss_motifRanking

    # Motif annotation (TFs)
    data(hg19_direct_motifAnnotation)
    direct_motifAnnotation <- hg19_direct_motifAnnotation
    data(hg19_inferred_motifAnnotation) # optional
    inferred_motifAnnotation <- hg19_inferred_motifAnnotation
  }
  allTFs <- unique(direct_motifAnnotation$allTFs, inferred_motifAnnotation$allTFs)

  ################################################################################
  # Step 1. Part 2: Creating TF modules (potential TF-targets)
  ################################################################################
  options(stringsAsFactors=FALSE) #TO DO: Better way to set it?

  if("1.2" %in% stepsToRun)
  {
    msg <- paste0(format(Sys.time(), "%H:%M"), "\tCreating TF modules")
    message(msg)

    # Convert the weight matrix into links:
    if(file.exists("int/1.5_GENIE3_linkList.RData"))
    {
      load("int/1.5_GENIE3_linkList.RData")
    }else{
      load("int/1.3_GENIE3_weightMatrix.RData")
      linkList <- getLinkList(weightMatrix, threshold=0.001) # (slighly faster)
      #linkList <- getLinkList(weightMatrix)
      colnames(linkList) <- c("TF", "Target", "weight")
      # order by weight
      linkList <- linkList[order(linkList[,"weight"], decreasing=TRUE),]
      save(linkList, file="int/1.5_GENIE3_linkList.RData")
    }


    dim(linkList)
    head(linkList)

    quantile(linkList$weight, probs=c(0.75, 0.90))
    plot(linkList$weight[1:1000000], type="l", ylim=c(0, max(linkList$weight)), main="Weight of the links",
         ylab="Weight", xlab="Links sorted decreasingly")
    abline(h=0.001, col="blue") # Threshold
    sum(linkList$weight>0.001)/nrow(linkList)

    linkList_001 <- linkList[which(linkList[,"weight"]>0.001),]
    # Number of links over the threshold:
    nrow(linkList_001)

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
    save(linkList_001_byTarget, file="int/1.5_linkList_001_byTarget.RData")

    nTopTfs <- c(5, 10, 50)
    nTopTfs <- setNames(nTopTfs, paste("top", nTopTfs, "perTarget", sep=""))

    library(reshape2); library(data.table)
    topTFsperTarget <- lapply(linkList_001_byTarget, function(llt) {
      nTFs <- nTopTfs[which(nTopTfs <= nrow(llt))]
      melt(lapply(nTFs, function(x) llt[1:x,"TF"]))
    })

    topTFsperTarget <- topTFsperTarget[which(!sapply(sapply(topTFsperTarget, nrow), is.null))]
    topTFsperTarget.asDf <-  data.frame(rbindlist(topTFsperTarget, idcol=TRUE))
    head(topTFsperTarget.asDf)
    colnames(topTFsperTarget.asDf) <- c("Target", "TF", "method")

    # Merge the all the gene-sets:
    tfModules.melted <- melt(tfModules)
    colnames(tfModules.melted) <- c("Target", "TF", "method")
    tfModules <- rbind(tfModules.melted, topTFsperTarget.asDf)

    save(tfModules, file="int/1.6_tfModules.RData")

    # Basic counts:
    rbind(nGeneSets=nrow(tfModules),
          nTFs=length(unique(tfModules$TF)),
          nTargets=length(unique(tfModules$Target)))

    ### Split into positive- and negative-correlated targets
    load("int/1.4_corrMat.RData")
    # Keep only correlation between TFs and potential targets
    tfs <- unique(tfModules$TF)
    corrMat <- corrMat[tfs,]

    # Split TF modules according to correlation
    tfModules_byTF <- split(tfModules, factor(tfModules$TF))
    tfModules_withCorr_byTF <- lapply(tfModules_byTF, function(tfGeneSets)
    {
      tf <- unique(tfGeneSets$TF)
      targets <- tfGeneSets$Target
      cbind(tfGeneSets, corr=c(as.numeric(corrMat[tf,targets] > 0.03) - as.numeric(corrMat[tf,targets] < -0.03)))
    })
    tfModules_withCorr <- data.frame(rbindlist(tfModules_withCorr_byTF))
    save(tfModules_withCorr, file="int/1.7_tfModules_withCorr.RData")
  }

  ################################################################################
  # Step 2. Identifying regulons (direct TF targets) based on DNA motif enrichment
  ################################################################################

  if("2" %in% stepsToRun)
  {
    if(!"1.2" %in% stepsToRun)  # If 1.t has not been run in this execution: Load its outputs
    {
      if(!file.exists("int/1.7_tfModules_withCorr.RData")){
        stop("The output from Step 1 should be stored in 'int/1.7_tfModules_withCorr.RData'")
      }
      load("int/1.7_tfModules_withCorr.RData")
    }

    msg <- paste0(format(Sys.time(), "%H:%M"), "\tStep 2. Identifying regulons")
    message(msg)

    ### Filter & format co-expression modules
    # Remove genes missing from RcisTarget databases
    #  (In case the input matrix wasn't already filtered)
    tfModules_withCorr <- tfModules_withCorr[which(as.character(tfModules_withCorr$TF) %in% allTFs),]
    geneInDb <- tfModules_withCorr$Target %in% motifRankings[["500bp"]]@rankings$rn
        missingGene <- sort(unique(tfModules_withCorr[which(!geneInDb),"Target"]))
        if(length(missingGene)>0) warning(paste0("Genes in co-expression modules not available in RcisTargetDatabases: ", paste(missingGene, collapse=", ")))
    tfModules_withCorr <- tfModules_withCorr[which(geneInDb),]

    # Targets with positive correlation
    tfModules_Selected <- tfModules_withCorr[which(tfModules_withCorr$corr==1),]

    # Add a column with the geneSet name (TF_method)
    tfModules_Selected <- cbind(tfModules_Selected, geneSetName=paste(tfModules_Selected$TF, tfModules_Selected$method, sep="_"))
    head(tfModules_Selected)

    # Split into tfModules (TF-modules, with several methods)
    tfModules <- split(tfModules_Selected$Target, tfModules_Selected$geneSetName)

    # Keep gene sets with at least 20 genes
    tfModules <- tfModules[which(lengths(tfModules)>=20)]

    # Add TF to the gene set (used in the following steps, careful if editing)
    tfModules <- setNames(lapply(names(tfModules), function(gsn) {
      tf <- strsplit(gsn, "_")[[1]][1]
      unique(c(tf, tfModules[[gsn]]))
    }), names(tfModules))
    save(tfModules, file="int/2.1_tfModules_forMotifEnrichmet.RData")

    tfModulesSummary <- t(sapply(strsplit(names(tfModules), "_"), function(x) x[1:2]))
    if(verbose) {
        message("tfModulesSummary:")
        print(sort(table(tfModulesSummary[,2])))
    }

    ################################################################
    ### 1. Calculate motif enrichment for each TF-module (Run RcisTarget)

    ### 1.1 Calculate enrichment
    msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Calculating AUC")
    message(msg)

    motifs_AUC <- lapply(motifRankings, function(ranking) calcAUC(tfModules, ranking, aucMaxRank=0.01*nrow(ranking), nCores=nCores, verbose=FALSE))
    save(motifs_AUC, file="int/2.2_motifs_AUC_500bp_10kbp.RData")
    # load("int/2.2_motifs_AUC_500bp_10kbp.RData")

    ### 1.2 Conver to table, filter by NES & add the TFs to which the motif is annotated
    # (For each database...)
    msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Adding motif annotation")
    message(msg)
    motifEnrichment <- lapply(motifs_AUC, function(aucOutput)
    {
      # Extract the TF of the gene-set name (i.e. MITF_w001):
      tf <- sapply(setNames(strsplit(rownames(aucOutput), "_"), rownames(aucOutput)), function(x) x[[1]])

      # Calculate NES and add motif annotation (provide tf in 'highlightTFs'):
      addMotifAnnotation(aucOutput, highlightTFs=tf, nesThreshold=3.0, digits=3,
                         motifAnnot_direct=direct_motifAnnotation,
                         motifAnnot_inferred=inferred_motifAnnotation)
    })

    # Merge both tables, adding a column that contains the 'motifDb'
    motifEnrichment <- do.call(rbind, lapply(names(motifEnrichment), function(dbName){
      cbind(motifDb=dbName, motifEnrichment[[dbName]])
    }))
    save(motifEnrichment, file="int/2.3_motifEnrichment.RData")
    msg <- paste0("Number of motifs in the initial enrichment: ", nrow(motifEnrichment))
    message(msg)

    ### 1.3 Keep only the motifs annotated to the initial TF
    motifEnrichment_selfMotifs <- motifEnrichment[which(motifEnrichment$TFinDB != ""),, drop=FALSE]
    save(motifEnrichment_selfMotifs, file="int/2.4_motifEnrichment_selfMotifs.RData")
    msg <- paste0("Number of motifs annotated to the initial TF: ", nrow(motifEnrichment_selfMotifs))
    message(msg)
    rm(motifEnrichment)

    ################################################################
    # 2. Prune targets
    # load("int/2.1_tfModules_forMotifEnrichmet.RData")
    # load("int/2.4_motifEnrichment_selfMotifs.RData")

    msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Prunning targets")
    message(msg)
    motifEnrichment_selfMotifs_wGenes <- lapply(names(motifRankings), function(motifDbName){
      addSignificantGenes(resultsTable=motifEnrichment_selfMotifs[motifDb==motifDbName],
                          geneSets=tfModules,
                          rankings=motifRankings[[motifDbName]],
                          maxRank=5000, method="aprox", nCores=nCores)
    })

    library(data.table)
    motifEnrichment_selfMotifs_wGenes <- rbindlist(motifEnrichment_selfMotifs_wGenes)
    save(motifEnrichment_selfMotifs_wGenes, file="int/2.5_motifEnrichment_selfMotifs_wGenes.RData")

    # Save as text:
    write.table(motifEnrichment_selfMotifs_wGenes, file="output/Step2_MotifEnrichment.tsv",
                sep="\t", quote=FALSE, row.names=FALSE)


    # load("int/2.5_motifEnrichment_selfMotifs_wGenes.RData")
    dim(motifEnrichment_selfMotifs_wGenes)
    motifEnrichment_selfMotifs_wGenes[order(NES,decreasing=TRUE)][1:5,-"enrichedGenes", with=F]

    ################################################################
    # Format regulons & save
    library(data.table)
    motifEnrichment.asIncidList <- apply(motifEnrichment_selfMotifs_wGenes, 1, function(oneMotifRow) {
      genes <- strsplit(oneMotifRow["enrichedGenes"], ";")[[1]]
      oneMotifRow <- data.frame(rbind(oneMotifRow), stringsAsFactors=FALSE)
      data.frame(oneMotifRow[rep(1, length(genes)),c("NES", "motif", "highlightedTFs", "TFinDB")], genes, stringsAsFactors = FALSE)
    })
    motifEnrichment.asIncidList <- rbindlist(motifEnrichment.asIncidList)
    colnames(motifEnrichment.asIncidList) <- c("NES", "motif", "TF", "annot", "gene")
    motifEnrichment.asIncidList <- data.frame(motifEnrichment.asIncidList, stringsAsFactors = FALSE)

    # Get targets for each TF, but keep info about best motif/enrichment
    # (directly annotated motifs are considered better)
    regulonTargetsInfo <- lapply(split(motifEnrichment.asIncidList, motifEnrichment.asIncidList$TF), function(tfTargets){
      # print(unique(tfTargets$TF))
      tfTable <- as.data.frame(do.call(rbind, lapply(split(tfTargets, tfTargets$gene), function(enrOneGene){
        directAnnot <- "**" %in% enrOneGene$annot
        enrOneGeneByAnnot <- enrOneGene
        if(directAnnot) enrOneGeneByAnnot <- enrOneGeneByAnnot[which(enrOneGene$annot == "**"),]
        bestMotif <- which.max(enrOneGeneByAnnot$NES)

        cbind(TF=unique(enrOneGene$TF), gene=unique(enrOneGene$gene), nMotifs=nrow(enrOneGene),
              bestMotif=as.character(enrOneGeneByAnnot[bestMotif,"motif"]), NES=as.numeric(enrOneGeneByAnnot[bestMotif,"NES"]),
              directAnnot=directAnnot)
      })), stringsAsFactors=FALSE)
      tfTable[order(tfTable$NES, decreasing = TRUE),]
    })
    rm(motifEnrichment.asIncidList)
    regulonTargetsInfo <- rbindlist(regulonTargetsInfo)
    colnames(regulonTargetsInfo) <- c("TF", "gene", "nMotifs", "bestMotif", "NES", "directAnnot")

    # Optional: Add Genie3 score
    load("int/1.5_GENIE3_linkList.RData")
    if("weight" %in% colnames(linkList))
    {
      linkList <- linkList[which(linkList$weight>=0.001),]  # TO DO: Will not work with GRNBOOST!
    }
    rownames(linkList) <- paste(linkList$TF, linkList$Target,sep="__")
    regulonTargetsInfo <- cbind(regulonTargetsInfo, Genie3Weight=linkList[paste(regulonTargetsInfo$TF, regulonTargetsInfo$gene,sep="__"),"weight"])

    save(regulonTargetsInfo, file="int/2.6_regulonTargetsInfo.RData")
    write.table(regulonTargetsInfo, file="output/Step2_regulonTargetsInfo.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    rm(linkList)

    # Split into regulons... (output: list TF --> targets)
    regulonTargetsInfo_splitByAnnot <- split(regulonTargetsInfo, regulonTargetsInfo$directAnnot)
    regulons <- sapply(split(regulonTargetsInfo_splitByAnnot[["TRUE"]], regulonTargetsInfo_splitByAnnot[["TRUE"]][,"TF"]), function(x) sort(as.character(unlist(x[,"gene"]))))
    regulons_extended <- sapply(split(regulonTargetsInfo_splitByAnnot[["FALSE"]],regulonTargetsInfo_splitByAnnot[["FALSE"]][,"TF"]), function(x) unname(x[,"gene"]))
    regulons_extended <- sapply(names(regulons_extended), function(tf) sort(unique(c(regulons[[tf]], regulons_extended[[tf]]))))
    names(regulons_extended) <- paste(names(regulons_extended), "_extended", sep="")
    regulons <- c(regulons, regulons_extended)
    save(regulons, file="int/2.6_regulons_asGeneSet.RData")


    # Number of regulons and summary of sizes:
    length(regulons)
    summary(lengths(regulons))

    # Save as incidence matrix (i.e. network)
    incidList <- melt(regulons)
    incidMat <- table(incidList[,2], incidList[,1])
    save(incidMat, file="int/2.6_regulons_asIncidMat.RData")
    rm(incidMat)
  }

  ################################################################################
  # Step 3. Analyzing the network activity in each individual cell
  ################################################################################

  if("3.1" %in% stepsToRun)
  {
    if(!"2" %in% stepsToRun)  # If 2 has not been run in this execution: Load its outputs TO DO: (something needed from 1?)
    {
      # load("int/2.2_motifs_AUC_500bp_10kbp.RData")
      load("int/2.6_regulons_asGeneSet.RData")
    #   if(!file.exists("int/")){
    #     # stop("The output from xxx should be stored in 'int/xxx.RData'")
    #   }
    }

    msg <- paste0(format(Sys.time(), "%H:%M"), "\tStep 3. Analyzing the network activity in each individual cell")
    message(msg)

    ## Prepare regulons
    regulons <- regulons[order(lengths(regulons), decreasing=TRUE)]
    regulons <- regulons[lengths(regulons)>=10]

    # Add the TF & rename
    regulons <- setNames(lapply(names(regulons), function(tf) sort(unique(c(gsub("_extended", "", tf), regulons[[tf]])))), names(regulons))
    names(regulons) <- paste(names(regulons), " (",lengths(regulons), "g)", sep="")

    save(regulons, file="int/3.0_regulons_forAUCell.RData")

    msg <- paste0("\nNumber of regulons to evaluate on cells: ", length(regulons),
                  "\nBiggest (non-extended) regulons: \n",
                  paste("\t", names(regulons)[which(!grepl("_extended", names(regulons)))][1:10], collapse="\n"))
    message(msg)



    ################################################################
    # AUCell
    library(AUCell)
    # 1. Create rankings
    pdf("output/Step3.1_aucellGenesStats.pdf")
    aucellRankings <- AUCell.buildRankings(exprMat, nCores=nCores, plotStats=TRUE)
    abline(v=aucellRankings@nGenesDetected["1%"], col="skyblue3", lwd=5, lty=3)
    dev.off()
    save(aucellRankings, file="int/3.1_aucellRankings.RData")

    # 2. Calculate AUC
    regulonAUC <- AUCell.calcAUC(regulons, aucellRankings, aucMaxRank=aucellRankings@nGenesDetected["1%"], nCores=nCores)

    # Order the modules by similarity, for easier exploration in the upcoming steps & save
    variableRegulons <- names(which(apply(getAuc(regulonAUC), 1, sd) > 0))
    reguDist <-as.dist(1-cor(t(getAuc(regulonAUC)[variableRegulons,]), method="spear"))
    reguClust <- hclust(reguDist, method="ward.D2")
    regulonClusters <- setNames(dynamicTreeCut::cutreeDynamic(reguClust, distM=as.matrix(reguDist), verbose = FALSE), reguClust$labels)
    regulonOrder <- reguClust$labels[reguClust$order]
    regulonOrder <- regulonOrder[order(regulonClusters[regulonOrder], decreasing = TRUE)]
    regulonAUC@matrix <- regulonAUC@matrix[regulonOrder,]
    save(regulonAUC, file="int/3.2_regulonAUC.RData")

    #### Overview of cell states according to module activity (tSNE on AUC)
    msg <- paste0(format(Sys.time(), "%H:%M"), "\tPlotting t-SNEs and histograms")
    message(msg)

    # Number of genes detected:
    nGenesPerCell <- apply(exprMat, 2, function(x) sum(x>0))
    colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
    cellColorNgenes <- setNames(adjustcolor(colorPal(10), alpha=.8)[as.numeric(cut(nGenesPerCell,breaks=10, right=F,include.lowest=T))], names(nGenesPerCell))
    save(cellColorNgenes, file="int/cellColorNgenes.RData")

    regulonAUC_subset <- subset(regulonAUC, onlyNonDirectExtended(rownames(regulonAUC)))

    # PCA-based t-SNE
    set.seed(seed)
    # (It is recommended to try different perplexity values)
    tsneAUC <- Rtsne::Rtsne(t(getAuc(regulonAUC_subset)), initial_dims=10, perplexity=10)
    rownames(tsneAUC$Y) <- colnames(regulonAUC_subset)
    colnames(tsneAUC$Y) <- c("tsne1", "tsne2")
    save(tsneAUC, file="int/3.3_tsneRegulonAUC_PCA.RData")

    if(showAlternativeTsnes)
    {
      # (More separation across clusters)
      tsneAUC <- Rtsne::Rtsne(t(getAuc(regulonAUC_subset)), initial_dims=50, perplexity=50)
      rownames(tsneAUC$Y) <- colnames(regulonAUC_subset)
      colnames(tsneAUC$Y) <- c("tsne1", "tsne2")
      save(tsneAUC, file="int/3.3_tsneRegulonAUC_PCA_50pc50perpl.RData")

      # Alternative: Distance-based t-SNE:
      corDist <- as.dist(1-cor(getAuc(regulonAUC_subset)))
      set.seed(seed)
      tsneAUC <- Rtsne::Rtsne(corDist, is_distance=TRUE, perplexity=10)
      rownames(tsneAUC$Y) <- labels(corDist)
      colnames(tsneAUC$Y) <- c("tsne1", "tsne2")
      save(tsneAUC, file="int/3.3_tsneRegulonAUC_Dist.RData")

      # Plot the preview of the t-SNEs (only colored by genes)
      pdf("output/Step3.1_AUCtSNEs.pdf")
      par(mfrow=c(2,2))
      load("int/3.3_tsneRegulonAUC_PCA.RData")
      tSNE <- tsneAUC$Y
      plot(tSNE, col=cellColorNgenes[rownames(tSNE)], pch=16, main="10 PCs, 10 perpl",
           xlab="t-SNE on regulon activity (AUC)", ylab="",
           sub=paste("Color: # detected genes"), axes=FALSE); box()

      load("int/3.3_tsneRegulonAUC_PCA_50pc50perpl.RData")
      tSNE <- tsneAUC$Y
      plot(tSNE, col=cellColorNgenes[rownames(tSNE)], pch=16, main="50 PCs, 50 perpl",
           xlab="t-SNE on regulon activity (AUC)", ylab="",
           sub=paste("Color: # detected genes"), axes=FALSE); box()

      load("int/3.3_tsneRegulonAUC_Dist.RData")
      tSNE <- tsneAUC$Y
      plot(tSNE, col=cellColorNgenes[rownames(tSNE)], pch=16, main="Correlation, 10 perpl",
           xlab="t-SNE on regulon activity (AUC)", ylab="",
           sub=paste("Color: # detected genes"), axes=FALSE); box()
      plot.new()

      for(varName in colnames(cellInfo))
      {
        cellColor <- setNames(colVars[[varName]][as.character(cellInfo[,varName])], rownames(cellInfo))

        load("int/3.3_tsneRegulonAUC_PCA.RData")
        tSNE <- tsneAUC$Y
        plot(tSNE, col=cellColor[rownames(tSNE)], pch=16, main="10 PCs, 10 perpl",
            xlab="t-SNE on regulon activity (AUC)", ylab="",
            sub=paste("Color:",varName), axes=FALSE); box()
        box()

        if(showAlternativeTsnes){
        load("int/3.3_tsneRegulonAUC_PCA_50pc50perpl.RData")
        tSNE <- tsneAUC$Y
        plot(tSNE, col=cellColor[rownames(tSNE)], pch=16, main="50 PCs, 50 perpl",
             xlab="t-SNE on regulon activity (AUC)", ylab="",
             sub=paste("Color:",varName), axes=FALSE); box()

        load("int/3.3_tsneRegulonAUC_Dist.RData")
        tSNE <- tsneAUC$Y
        plot(tSNE, col=cellColor[rownames(tSNE)], pch=16, main="Correlation, 10 perpl",
             xlab="t-SNE on regulon activity (AUC)", ylab="",
             sub=paste("Color:",varName), axes=FALSE); box()
        }
      }

      dev.off()
    }

    ################################################################
    # Plot AUC histograms

    load("int/3.3_tsneRegulonAUC_PCA.RData")
    tSNE <- tsneAUC$Y
    par(mfrow=c(1,2))

    Cairo::CairoPDF("output/Step3.1_RegulonActivity_AUCtSNE.pdf", width=20, height=5)
    par(mfrow=c(1,4))

    # tSNE (colored by number of genes detected per cell)
    plot(tSNE, col=cellColorNgenes[rownames(tSNE)], pch=16, main="nGenes", sub="t-SNE on regulon activity (AUC)")

    # Other known properties:
    nVars <- 0
    if(!is.null(cellInfo))
    {
      nVars <- ncol(cellInfo)

      for(varName in colnames(cellInfo))
      {
        if(is.null(colVars[[varName]])) {
          colVars[[varName]] <- setNames(rainbow(length(unique(cellInfo[,varName]))), unique(cellInfo[,varName]))
          cellColor <- setNames(colVars[[varName]][cellInfo[,varName]], rownames(cellInfo))
        }
        if(!is.null(colVars[[varName]]))
          cellColor <- setNames(colVars[[varName]][as.character(cellInfo[,varName])],
                                rownames(cellInfo))

        plot(tSNE, col=cellColor[rownames(tSNE)], pch=16, main=varName, sub="t-SNE on regulon activity (AUC)")
      }
    }

    for( i in seq_len(4 - (nVars+1 %% 4))) # fill remaining slots in page
    {
      plot.new()
    }


    # Plot module activity, thresholds & assignment:
    cells_AUCellThresholds <- plot_aucTsne(tSNE=tSNE, exprMat=exprMat, regulonAUC=regulonAUC, alphaOff=0.1)
    dev.off()
    save(cells_AUCellThresholds, file="int/3.4_AUCellThresholds.RData")

    ################################################################
    # Binary regulon activity matrix (Active regulons per cell)

    # Get cells assigned to each regulon
    regulonsCells <- lapply(cells_AUCellThresholds, function(x) x$assignment)

    ### Save threshold info as text (e.g. to edit/modify...)
    trhAssignment <- sapply(cells_AUCellThresholds, function(x) unname(x$aucThr$selected))
    commentsThresholds <- sapply(cells_AUCellThresholds, function(x) unname(x$aucThr$comment))

    table2edit <- cbind(regulon=names(trhAssignment),
                        threshold=trhAssignment,
                        nCellsAssigned=lengths(regulonsCells)[names(trhAssignment)],
                        AUCellComment=commentsThresholds,
                        nGenes=gsub("[\\(g\\)]", "", regmatches(names(cells_AUCellThresholds), gregexpr("\\(.*?\\)", names(cells_AUCellThresholds)))),
                        clusteringOrder=1:length(trhAssignment),
                        clusterGroup=regulonClusters[names(trhAssignment)],
                        onlyNonDirectExtended=(names(trhAssignment) %in% onlyNonDirectExtended(names(trhAssignment))),
                        personalNotes="")
    write.table(table2edit, file="int/3.5_1_AUCellThresholds.txt", row.names=F, quote=F, sep="\t")
    rm(trhAssignment)
  }

  if("3.2" %in% stepsToRun)
  {
    if(!"3.1" %in% stepsToRun)  # If 3.1 has not been run in this execution: Load its outputs TO DO: (something needed from 1 or 2?)
    {
      #   if(!file.exists("int/")){
      #     # stop("The output from xxx should be stored in 'int/xxx.RData'")
      #   }
      # regulonsCells
    }

    ### Convert to matrix (regulons with zero assigned cells are lost)
    regulonActivity <- reshape2::melt(regulonsCells)
    binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
    class(binaryRegulonActivity) <- "matrix"
    save(binaryRegulonActivity, file="int/3.6_binaryRegulonActivity.RData")

    # Keep only non-duplicated thresholds
    # (e.g. only "extended" regulons if there is not a regulon based on direct annotation)
    binaryRegulonActivity_nonDupl <- binaryRegulonActivity[which(rownames(binaryRegulonActivity) %in% onlyNonDirectExtended(rownames(binaryRegulonActivity))),]
    save(binaryRegulonActivity_nonDupl, file="int/3.7_binaryRegulonActivity_nonDupl.RData")

    minCells <- ncol(exprMat) * .01
    msg <- paste0("Binary regulon activity: ",
                  nrow(binaryRegulonActivity_nonDupl), " TF regulons x ",
                  ncol(binaryRegulonActivity), " cells.\n(",
                  nrow(binaryRegulonActivity), " regulons including 'extended' versions)\n",
                  sum(rowSums(binaryRegulonActivity_nonDupl)>minCells),
                  " regulons are active in more than 1% (", minCells, ") cells.")
    message(msg)


    ### Boxplot active cells per regulon
    Cairo::CairoPDF("output/Step3.2_BoxplotActiveCellsRegulon.pdf", width=6, height=6)
    par(mfrow=c(1,2))
    boxplot(rowSums(binaryRegulonActivity_nonDupl), main="nCells per regulon",
            sub='number of cells \nthat have the regulon active',
            col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)
    boxplot(colSums(binaryRegulonActivity_nonDupl), main="nRegulons per Cell",
            sub='number of regulons \nactive per cell',
            col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)
    dev.off()

    ################################################################################
    # Binary activity heatmap

    ### Select regulons:
    regulonSelection <- list()

    # All regulons.
    regulonSelection[["All regulons \n (including duplicated regulons)"]] <- rownames(binaryRegulonActivity)

    # Active in > 1% cells
    regMinCells <- names(which(rowSums(binaryRegulonActivity_nonDupl) > minCells))
    regulonSelection[["Regulons active in more than 1% of cells"]] <- regMinCells

    # Correlation across regulons (based on binary cell activity)
    reguCor <- cor(t(binaryRegulonActivity_nonDupl[regMinCells,]))
    diag(reguCor) <- 0

    # Regulons that co-ocurr in similar cells. If a regulon is relevant by itself it will not be shown, also check the regulons ignored.
    corrRegs <- names(which(rowSums(abs(reguCor) > 0.30) > 0))
    regulonSelection[["Regulons with any other regulon correlated\n with abs(cor)>0.30 \n(and active in at least 1% of cells)"]]  <- corrRegs
    missingRegs <- rownames(binaryRegulonActivity_nonDupl)[which(!rownames(binaryRegulonActivity_nonDupl) %in% corrRegs)]
    regulonSelection[["Regulons no other regulons correlated\n with abs(cor)>0.30 \n or active in fewer than 1% of cells"]]  <- missingRegs
    save(regulonSelection,file="int/3.8_regulonSelections.RData")

    ## Set regulon order (only plotting most correlated regulons)
    binaryRegulonOrder <- hclust(as.dist(1-reguCor[corrRegs,corrRegs]))
    binaryRegulonOrder <- binaryRegulonOrder$labels[binaryRegulonOrder$order]
    save(binaryRegulonOrder,file="int/3.9_binaryRegulonOrder.RData")

    ### Plot heatmap:
    for(i in seq_len(length(regulonSelection)))
    {
      selRegs <- names(regulonSelection)[i]
      binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]],,drop=FALSE]
      NMF::aheatmap(binaryMat, scale="none", revC=TRUE, main=selRegs,
                    annCol=cellInfo[colnames(binaryMat),, drop=FALSE],
                    annColor=colVars,
                    color = c("white", "black"),
                    filename=paste0("output/Step3.3_binaryRegulonActivity_Heatmap_",i,".pdf"))
    }
  }
  ################################################################################
  # Step 4. Clustering. Here: Preview of t-SNEs
  ################################################################################
  if("4" %in% stepsToRun)
  {
    if(!"3.2" %in% stepsToRun)  # If 3.2 has not been run in this execution: Load its outputs TO DO: (something needed from previous?)
    {
      #   if(!file.exists("int/")){
      #     # stop("The output from xxx should be stored in 'int/xxx.RData'")
      #   }
      load("int/cellColorNgenes.RData")
      load("int/3.7_binaryRegulonActivity_nonDupl.RData")
      load("int/3.9_binaryRegulonOrder.RData")
    }

    library(Rtsne)
    tBinaryAct <- t(binaryRegulonActivity_nonDupl)

    ##################################
    # PCA based t-SNE
    set.seed(123)
    tBinaryAct_jitter <- jitter(tBinaryAct, factor=1)
    tsneBinaryActivity_PCA <- Rtsne(tBinaryAct_jitter, initial_dims=5, perplexity=30)
    rownames(tsneBinaryActivity_PCA$Y) <- rownames(tBinaryAct_jitter)
    colnames(tsneBinaryActivity_PCA$Y) <- c("tsne2", "tsne1")
    tsneBinaryActivity_PCA$Y <- tsneBinaryActivity_PCA$Y[,c("tsne1", "tsne2")]
    save(tsneBinaryActivity_PCA, file="int/4.1_tsneBinaryActivity_5PC.RData")

    ##################################
    # PCA based t-SNE
    set.seed(123)
    tBinaryAct_jitter <- jitter(tBinaryAct, factor=1)
    tsneBinaryActivity_PCA <- Rtsne(tBinaryAct_jitter, initial_dims=50, perplexity=30)
    rownames(tsneBinaryActivity_PCA$Y) <- rownames(tBinaryAct_jitter)
    colnames(tsneBinaryActivity_PCA$Y) <- c("tsne2", "tsne1")
    tsneBinaryActivity_PCA$Y <- tsneBinaryActivity_PCA$Y[,c("tsne1", "tsne2")]
    save(tsneBinaryActivity_PCA, file="int/4.1_tsneBinaryActivity_50PC.RData")

    ##################################
    # Distance-based t-SNE
    corDist <- as.dist(1-cor(t(tBinaryAct)))
    set.seed(123)
    tsneBinaryActivity_Dist <- Rtsne(corDist, is_distance=TRUE, perplexity=30)
    rownames(tsneBinaryActivity_Dist$Y) <- labels(corDist)
    colnames(tsneBinaryActivity_Dist$Y) <- c("tsne1", "tsne2")
    save(tsneBinaryActivity_Dist, file="int/4.1_tsneBinaryActivity_Dist.RData")

    ##################################
    # Plot the different versions of the t-SNEs... to compare/decide what is best
    tSNEs_binary <- list()
    tSNEs_binary[["Dist"]] <- tsneBinaryActivity_Dist$Y
    load("int/4.1_tsneBinaryActivity_5PC.RData")
    tSNEs_binary[["5PC"]] <- tsneBinaryActivity_PCA$Y
    load("int/4.1_tsneBinaryActivity_50PC.RData")
    tSNEs_binary[["50PC"]] <- tsneBinaryActivity_PCA$Y

    for(tsneName in names(tSNEs_binary))
    {
      tSNE_binary <- tSNEs_binary[[tsneName]]
      # to detect most likely stable states (higher-density areas in the t-SNE)
      # We can also generate a landscape out of the tSNE data, in which zones with higher density can be spotted. These zones would correspond to 'stable states':
      library(KernSmooth)
      library(RColorBrewer)
      dens2d <- bkde2D(tSNE_binary, 1)$fhat

      Cairo::CairoPDF(paste0("output/Step4.1_tsneModuleActivity_",tsneName,".pdf"), width=15, height=5)
      par(mfrow=c(1,3))
      # nGenes
      plot(tSNE_binary, col=cellColorNgenes[rownames(tSNE_binary)], pch=16)
      # density
      image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
      contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

      # Known phenotype:
      if(!is.null(cellInfo))
      {
        nVars <- ncol(cellInfo)
        for(varName in colnames(cellInfo))
        {
          cellColor <- setNames(colVars[[varName]][as.character(cellInfo[,varName])], rownames(cellInfo))
          plot(tSNE_binary, col=cellColor[rownames(tSNE_binary)], pch=16, main=varName, sub="t-SNE on Binary regulon activity",
               xlab="", ylab="",axes = FALSE)
        }
      }
      # legend(10, 25, names(colVars[[varName]]), fill=colVars[[varName]], cex=.7, bty="n")
      for(i in seq_len(3 - ((nVars+2) %% 3))) # fill remaining slots in page
      {
        plot.new()
      }
      dev.off()
    }

    ##################################
    ## Visualize module activity on one of these t-SNEs
    ## Here we are plotting it on 50pcs (it could also be done on any other tSNE)
    load("int/4.1_tsneBinaryActivity_50PC.RData")
    tSNE_binary <- tsneBinaryActivity_PCA$Y

    Cairo::CairoPDF("output/Step4.2_tsneBinaryActivity_50PC_BinaryRegulons.pdf", width=20, height=15)
    par(mfrow=c(4,6))
    cells_trhAssignment <- plot_aucTsne(tSNE=tSNE_binary, exprMat=exprMat, regulonAUC=tBinaryAct[binaryRegulonOrder,], cex=1.5, plots="binary", thresholds=0)
    dev.off()

    Cairo::CairoPDF("output/Step4.2_tsneBinaryActivity_50PC_AUCRegulons.pdf", width=20, height=15)
    par(mfrow=c(4,6))
    cells_trhAssignment <- plot_aucTsne(tSNE=tSNE_binary, exprMat=exprMat, regulonAUC=regulonAUC[binaryRegulonOrder,], cex=1.5, plots="AUC", thresholds=0)
    dev.off()


    Cairo::CairoPDF("output/Step4.2_tsneBinaryActivity_50PC_allPlots.pdf", width=20, height=5)
    par(mfrow=c(1,4))
    cells_trhAssignment <- plot_aucTsne(tSNE=tSNE_binary, exprMat=exprMat, regulonAUC=regulonAUC[binaryRegulonOrder,],
                                        alphaOff=0.1, thresholds=cells_AUCellThresholds[binaryRegulonOrder])
    dev.off()
  }
  sink()
}
