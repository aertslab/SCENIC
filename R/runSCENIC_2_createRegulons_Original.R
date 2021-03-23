#' ################################################################################
#' # Step 2. Identifying regulons (direct TF targets) based on DNA motif enrichment
#' ################################################################################
#' 
#' 
#' #' @title runSCENIC_2_createRegulons
#' #' @description Step 2: RcisTarget (prune co-expression modules using TF-motif enrichment analysis)
#' #' @param scenicOptions Fields used: TODO
#' #' @param minGenes Minimum size of co-expression gene set (default: 20 genes)
#' #' @param coexMethod Allows to select the method(s) used to generate the co-expression modules
#' #' @return The output is written in the folders 'int' and 'ouput'
#' #' @details See the detailed vignette explaining the internal steps.
#' #' @examples 
#' #' scenicOptions <- readRDS("int/scenicOptions.Rds")
#' #' # In case any settings need to be modified:
#' #' scenicOptions@settings$nCores <- 20
#' #' scenicOptions@inputDatasetInfo$org <- "mgi" 
#' #' 
#' #' runSCENIC_2_createRegulons(scenicOptions)
#' #' @export
#' runSCENIC_2_createRegulons <- function(scenicOptions, minGenes=20, coexMethod=NULL)
#' {
#'   nCores <- getSettings(scenicOptions, "nCores")
#'   tfModules_asDF <- loadInt(scenicOptions, "tfModules_asDF")
#'   if(!is.null(coexMethod)) tfModules_asDF <- tfModules_asDF[which(tfModules_asDF$method %in% coexMethod),]
#'   if(nrow(tfModules_asDF)==0) stop("The co-expression modules are empty.")
#'   
#'   # Set cores for RcisTarget::addMotifAnnotation(). The other functions use foreach package.
#'   if("BiocParallel" %in% installed.packages()) library(BiocParallel); register(MulticoreParam(nCores), default=TRUE) 
#'   
#'   msg <- paste0(format(Sys.time(), "%H:%M"), "\tStep 2. Identifying regulons")
#'   if(getSettings(scenicOptions, "verbose")) message(msg)
#' 
#'   ### Check org and load DBs
#'   if(is.na(getDatasetInfo(scenicOptions, "org"))) stop('Please provide an organism (scenicOptions@inputDatasetInfo$org).')
#'   library(AUCell)
#'   library(RcisTarget)
#'   motifAnnot <- getDbAnnotations(scenicOptions)
#'   
#'   if(is.null(names(getSettings(scenicOptions, "dbs")))) 
#'   {
#'     names(scenicOptions@settings$"dbs") <- scenicOptions@settings$"dbs"
#'     tmp <- sapply(strsplit(getSettings(scenicOptions, "dbs"),"-", fixed=T), function(x) x[grep("bp|kb",x)])
#'     if(all(lengths(tmp)>0)) names(scenicOptions@settings$"dbs") <- tmp
#'   }
#'   
#'   loadAttempt <- sapply(getDatabases(scenicOptions), dbLoadingAttempt)
#'   if(any(!loadAttempt)) stop("It is not possible to load the following databses: \n",
#'                                 paste(dbs[which(!loadAttempt)], collapse="\n"))
#'   
#'   genesInDb <- unique(unlist(lapply(getDatabases(scenicOptions), function(x)
#'     names(feather::feather_metadata(x)[["types"]]))))
#'   
#'   ### Filter & format co-expression modules
#'   # Remove genes missing from RcisTarget databases
#'   #  (In case the input matrix wasn't already filtered)
#'   tfModules_asDF$TF <- as.character(tfModules_asDF$TF)
#'   tfModules_asDF$Target <- as.character(tfModules_asDF$Target)
#'   allTFs <- getDbTfs(scenicOptions)
#'   tfModules_asDF <- tfModules_asDF[which(tfModules_asDF$TF %in% allTFs),]
#'   geneInDb <- tfModules_asDF$Target %in% genesInDb
#'       missingGene <- sort(unique(tfModules_asDF[which(!geneInDb),"Target"]))
#'       if(length(missingGene)>0) 
#'         warning(paste0("Genes in co-expression modules not available in RcisTargetDatabases: ", 
#'                                                paste(missingGene, collapse=", ")))
#'   tfModules_asDF <- tfModules_asDF[which(geneInDb),]
#' 
#'   # Targets with positive correlation
#'   tfModules_Selected <- tfModules_asDF[which(tfModules_asDF$corr==1),]
#' 
#'   # Add a column with the geneSet name (TF_method)
#'   tfModules_Selected <- cbind(tfModules_Selected, geneSetName=paste(tfModules_Selected$TF, tfModules_Selected$method, sep="_"))
#'   tfModules_Selected$geneSetName <- factor(as.character(tfModules_Selected$geneSetName))
#'   # head(tfModules_Selected)
#'   allGenes <- unique(tfModules_Selected$Target)
#' 
#'   # Split into tfModules (TF-modules, with several methods)
#'   tfModules <- split(tfModules_Selected$Target, tfModules_Selected$geneSetName)
#' 
#'   # Add TF to the gene set (used in the following steps, careful if editing)
#'   tfModules <- setNames(lapply(names(tfModules), function(gsn) {
#'     tf <- strsplit(gsn, "_")[[1]][1]
#'     unique(c(tf, tfModules[[gsn]]))
#'   }), names(tfModules))
#'   
#'   # Keep gene sets with at least 'minGenes' genes
#'   tfModules <- tfModules[which(lengths(tfModules)>=minGenes)]
#'   saveRDS(tfModules, file=getIntName(scenicOptions, "tfModules_forEnrichment")) #TODO as geneset? & previous step?
#' 
#'   if(getSettings(scenicOptions, "verbose")) {
#'       tfModulesSummary <- t(sapply(strsplit(names(tfModules), "_"), function(x) x[1:2]))
#'       message("tfModulesSummary:")
#'       print(sort(table(tfModulesSummary[,2])))
#'   }
#' 
#'   ################################################################
#'   ### 1. Calculate motif enrichment for each TF-module (Run RcisTarget)
#' 
#'   ### 1.1 Calculate enrichment
#'   msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Calculating AUC")
#'   if(getSettings(scenicOptions, "verbose")) message(msg)
#' 
#'   motifs_AUC <- lapply(getDatabases(scenicOptions), function(rnkName) {
#'     ranking <- importRankings(rnkName, columns=allGenes)
#'     message("Scoring database: ", ranking@description)
#'     RcisTarget::calcAUC(tfModules, ranking, aucMaxRank=0.03*getNumColsInDB(ranking), nCores=nCores, verbose=FALSE)})
#'   saveRDS(motifs_AUC, file=getIntName(scenicOptions, "motifs_AUC"))
#'   
#'   ### 1.2 Convert to table, filter by NES & add the TFs to which the motif is annotated
#'   # (For each database...)
#'   msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Adding motif annotation")
#'   message(msg)
#'   motifEnrichment <- lapply(motifs_AUC, function(aucOutput)
#'   {
#'     # Extract the TF of the gene-set name (i.e. MITF_w001):
#'     tf <- sapply(setNames(strsplit(rownames(aucOutput), "_"), rownames(aucOutput)), function(x) x[[1]])
#'     
#'     # Calculate NES and add motif annotation (provide tf in 'highlightTFs'):
#'     addMotifAnnotation(aucOutput, 
#'                        nesThreshold=3, digits=3, 
#'                        motifAnnot=motifAnnot,
#'                        motifAnnot_highConfCat=c("directAnnotation", "inferredBy_Orthology"),
#'                        motifAnnot_lowConfCat=c("inferredBy_MotifSimilarity",
#'                                                  "inferredBy_MotifSimilarity_n_Orthology"), 
#'                        highlightTFs=tf)
#'   })
#' 
#'   # Merge both tables, adding a column that contains the 'motifDb'
#'   motifEnrichment <- do.call(rbind, lapply(names(motifEnrichment), function(dbName){
#'     cbind(motifDb=dbName, motifEnrichment[[dbName]])
#'   }))
#'   saveRDS(motifEnrichment, file=getIntName(scenicOptions, "motifEnrichment_full"))
#'   msg <- paste0("Number of motifs in the initial enrichment: ", nrow(motifEnrichment))
#'   if(getSettings(scenicOptions, "verbose")) message(msg)
#' 
#'   ### 1.3 Keep only the motifs annotated to the initial TF
#'   motifEnrichment_selfMotifs <- motifEnrichment[which(motifEnrichment$TFinDB != ""),, drop=FALSE]
#'   msg <- paste0("Number of motifs annotated to the matching TF: ", nrow(motifEnrichment_selfMotifs))
#'   if(getSettings(scenicOptions, "verbose")) message(msg)
#'   rm(motifEnrichment)
#' 
#'   if(nrow(motifEnrichment_selfMotifs)==0) 
#'     stop("None of the co-expression modules present enrichment of the TF motif: There are no regulons.")
#'   
#'   ################################################################
#'   # 2. Prune targets
#'   msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Pruning targets")
#'   if(getSettings(scenicOptions, "verbose")) message(msg)
#'     
#'   dbNames <- getDatabases(scenicOptions)
#'   motifEnrichment_selfMotifs_wGenes <- lapply(names(dbNames), function(motifDbName){
#'     ranking <- importRankings(dbNames[motifDbName], columns=allGenes)
#'     addSignificantGenes(resultsTable=motifEnrichment_selfMotifs[motifEnrichment_selfMotifs$motifDb==motifDbName,],
#'                         geneSets=tfModules,
#'                         rankings=ranking,
#'                         maxRank=5000, method="aprox", nCores=nCores)
#'   })
#'   
#'   suppressPackageStartupMessages(library(data.table))
#'   motifEnrichment_selfMotifs_wGenes <- rbindlist(motifEnrichment_selfMotifs_wGenes)
#'   saveRDS(motifEnrichment_selfMotifs_wGenes, file=getIntName(scenicOptions, "motifEnrichment_selfMotifs_wGenes"))
#'   
#'   if(getSettings(scenicOptions, "verbose")) 
#'   {
#'     # TODO messages/print
#'     message(format(Sys.time(), "%H:%M"), "\tNumber of motifs that support the regulons: ", nrow(motifEnrichment_selfMotifs_wGenes))
#'     motifEnrichment_selfMotifs_wGenes[order(motifEnrichment_selfMotifs_wGenes$NES,decreasing=TRUE),][1:5,(1:ncol(motifEnrichment_selfMotifs_wGenes)-1), with=F] 
#'   }
#'   
#'   # Save as text:
#'   if(!file.exists("output")) dir.create("output") 
#'   write.table(motifEnrichment_selfMotifs_wGenes, file=getOutName(scenicOptions, "s2_motifEnrichment"),
#'               sep="\t", quote=FALSE, row.names=FALSE)
#'   if("DT" %in% installed.packages() && nrow(motifEnrichment_selfMotifs_wGenes)>0)
#'   {
#'     nvm <- tryCatch({
#'       colsToShow <- c("motifDb", "logo", "NES", "geneSet", "TF_highConf", "TF_lowConf")
#'       motifEnrichment_2html <- viewMotifs(motifEnrichment_selfMotifs_wGenes, colsToShow=colsToShow, options=list(pageLength=100))
#'       
#'       fileName <- getOutName(scenicOptions, "s2_motifEnrichmentHtml")
#'       
#'       dirName <- dirname(fileName)
#'       fileName <- basename(fileName)
#'       suppressWarnings(DT::saveWidget(motifEnrichment_2html, fileName))
#'       file.rename(fileName, file.path(dirName, fileName))
#'       if(getSettings(scenicOptions, "verbose")) message("\tPreview of motif enrichment saved as: ", file.path(dirName, fileName))
#'     }, error = function(e) print(e$message))
#'   }
#' 
#'   ################################################################
#'   # Format regulons & save
#'   motifEnrichment.asIncidList <- apply(motifEnrichment_selfMotifs_wGenes, 1, function(oneMotifRow) {
#'     genes <- strsplit(oneMotifRow["enrichedGenes"], ";")[[1]]
#'     oneMotifRow <- data.frame(rbind(oneMotifRow), stringsAsFactors=FALSE)
#'     data.frame(oneMotifRow[rep(1, length(genes)),c("NES", "motif", "highlightedTFs", "TFinDB")], genes, stringsAsFactors = FALSE)
#'   })
#'   motifEnrichment.asIncidList <- rbindlist(motifEnrichment.asIncidList)
#'   colnames(motifEnrichment.asIncidList) <- c("NES", "motif", "TF", "annot", "gene")
#'   motifEnrichment.asIncidList <- data.frame(motifEnrichment.asIncidList, stringsAsFactors = FALSE)
#' 
#'   # Get targets for each TF, but keep info about best motif/enrichment
#'   # (directly annotated motifs are considered better)
#'   regulonTargetsInfo <- lapply(split(motifEnrichment.asIncidList, motifEnrichment.asIncidList$TF), function(tfTargets){
#'     # print(unique(tfTargets$TF))
#'     tfTable <- as.data.frame(do.call(rbind, lapply(split(tfTargets, tfTargets$gene), function(enrOneGene){
#'       highConfAnnot <- "**" %in% enrOneGene$annot
#'       enrOneGeneByAnnot <- enrOneGene
#'       if(highConfAnnot) enrOneGeneByAnnot <- enrOneGeneByAnnot[which(enrOneGene$annot == "**"),]
#'       bestMotif <- which.max(enrOneGeneByAnnot$NES)
#' 
#'       cbind(TF=unique(enrOneGene$TF), gene=unique(enrOneGene$gene), nMotifs=nrow(enrOneGene),
#'             bestMotif=as.character(enrOneGeneByAnnot[bestMotif,"motif"]), NES=as.numeric(enrOneGeneByAnnot[bestMotif,"NES"]),
#'             highConfAnnot=highConfAnnot)
#'     })), stringsAsFactors=FALSE)
#'     tfTable[order(tfTable$NES, decreasing = TRUE),]
#'   })
#'   rm(motifEnrichment.asIncidList)
#'   regulonTargetsInfo <- rbindlist(regulonTargetsInfo)
#'   colnames(regulonTargetsInfo) <- c("TF", "gene", "nMotifs", "bestMotif", "NES", "highConfAnnot")
#' 
#'   # Optional: Add Genie3 score
#'   linkList <- loadInt(scenicOptions, "genie3ll", ifNotExists="null")
#'   if(!is.null(linkList) & ("weight" %in% colnames(linkList)))
#'   {
#'     if(is.data.table(linkList)) linkList <- as.data.frame(linkList)
#'     
#'     uniquePairs <- nrow(unique(linkList[,c("TF", "Target")]))
#'     if(uniquePairs == nrow(linkList)) {
#'       linkList <- linkList[which(linkList$weight>=getSettings(scenicOptions, "modules/weightThreshold")),]  # TODO: Will not work with GRNBOOST!
#'       rownames(linkList) <- paste(linkList$TF, linkList$Target,sep="__")
#'       regulonTargetsInfo <- cbind(regulonTargetsInfo, Genie3Weight=linkList[paste(regulonTargetsInfo$TF, regulonTargetsInfo$gene,sep="__"),"weight"])
#'     }else {
#'       warning("There are duplicated regulator-target (gene id/name) pairs in the co-expression link list.",
#'               "\nThe co-expression weight was not added to the regulonTargetsInfo table.")
#'     }
#'   }else warning("It was not possible to add the weight to the regulonTargetsInfo table.")
#' 
#'   saveRDS(regulonTargetsInfo, file=getIntName(scenicOptions, "regulonTargetsInfo"))
#'   
#'   write.table(regulonTargetsInfo, file=getOutName(scenicOptions, "s2_regulonTargetsInfo"),
#'               sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
#'   rm(linkList)
#' 
#'   # Split into regulons... (output: list TF --> targets)
#'   regulonTargetsInfo_splitByAnnot <- split(regulonTargetsInfo, regulonTargetsInfo$highConfAnnot)
#'   regulons <- NULL
#'   if(!is.null(regulonTargetsInfo_splitByAnnot[["TRUE"]]))
#'   {
#'     regulons <- lapply(split(regulonTargetsInfo_splitByAnnot[["TRUE"]], regulonTargetsInfo_splitByAnnot[["TRUE"]][,"TF"]), function(x) sort(as.character(unlist(x[,"gene"]))))
#'   }
#'   regulons_extended <- NULL
#'   if(!is.null(regulonTargetsInfo_splitByAnnot[["FALSE"]]))
#'   {
#'     regulons_extended <- lapply(split(regulonTargetsInfo_splitByAnnot[["FALSE"]],regulonTargetsInfo_splitByAnnot[["FALSE"]][,"TF"]), function(x) unname(unlist(x[,"gene"])))
#'     regulons_extended <- setNames(lapply(names(regulons_extended), function(tf) sort(unique(c(regulons[[tf]], unlist(regulons_extended[[tf]]))))), names(regulons_extended))
#'     names(regulons_extended) <- paste(names(regulons_extended), "_extended", sep="")
#'   }
#'   regulons <- c(regulons, regulons_extended)
#'   saveRDS(regulons, file=getIntName(scenicOptions, "regulons"))
#'   
#'   # Save as incidence matrix (i.e. network)
#'   incidList <- reshape2::melt(regulons)
#'   incidMat <- table(incidList[,2], incidList[,1])
#'   saveRDS(incidMat, file=getIntName(scenicOptions, "regulons_incidMat"))
#'   rm(incidMat)
#'   #TODO NMF::aheatmap(incidMat)
#'   
#'   if(getSettings(scenicOptions, "verbose")) 
#'   {
#'     # Number of regulons and summary of sizes:
#'     length(regulons) # TODO
#'     summary(lengths(regulons))
#'   }
#' }
#' 
#' 
#' #' @title getDbAnnotations
#' #' @description Loads the motif annotation
#' #' @param scenicOptions Fields used: 
#' #' If scenicOptions@settings$db_annotFiles is set, it will load these files.
#' #' Otherwise, will load the default RcisTarget annotations based on 'scenicOptions@inputDatasetInfo$org', and 'scenicOptions@settings$db_mcVersion'
#' #' @return The motif annotations
#' #' @examples 
#' #' getDbAnnotations(scenicOptions)
#' #' @export 
#' getDbAnnotations <- function(scenicOptions)
#' {
#'   dbAnnotFiles <- scenicOptions@settings$db_annotFiles
#'   if(!is.null(dbAnnotFiles))
#'   {
#'     motifAnnotations <- NULL
#'     for(annotPath in dbAnnotFiles)
#'     {
#'       motifAnnot <- data.table::fread(annotPath) #; head(motifAnnot)
#'       motifAnnot$annotationSource <- factor(motifAnnot$annotationSource)
#'       colnames(motifAnnot)[1]<-"motif" # for now...
#'       motifAnnotations <- rbind(motifAnnotations, motifAnnot)
#'     }
#'   } else { # Default RcisTarget annotations
#'     if(is.na(getDatasetInfo(scenicOptions, "org"))) stop('Please provide an organism (scenicOptions@inputDatasetInfo$org).')
#'     org <- getDatasetInfo(scenicOptions, "org")
#'     if(is.na(org)) stop("Please provide an organism (scenicOptions@inputDatasetInfo$org).")
#'     if(!org %in% c("hgnc", "mgi", "dmel")) stop("Organism not recognized (scenicOptions@inputDatasetInfo$org).")
#'     
#'     if(org=="hgnc") motifAnnotName <- "motifAnnotations_hgnc"
#'     if(org=="mgi") motifAnnotName <- "motifAnnotations_mgi"
#'     if(org=="dmel") motifAnnotName <- "motifAnnotations_dmel"
#'     
#'     if(!is.null(scenicOptions@settings$db_mcVersion))
#'     {
#'       if(scenicOptions@settings$db_mcVersion=="v8") motifAnnotName <- paste0(motifAnnotName, "_v8")
#'     }
#'     
#'     library(RcisTarget) # Lazyload
#'     #data(package="RcisTarget", verbose = T)
#'     data(list=motifAnnotName, package="RcisTarget", verbose = FALSE)
#'     motifAnnotations <- eval(as.name(motifAnnotName))
#'   }
#'   
#'   return(motifAnnotations)
#' }
#' 
#' #' @title getDbTfs
#' #' @description Provides the list of transcription factors in the RcisTarget databases for the given organism
#' #' @param scenicOptions Fields used: 'scenicOptions@inputDatasetInfo$org'
#' #' @return List of transcription factors in the databases.
#' #' @examples 
#' #' getDbTfs(scenicOptions)
#' #' @export 
#' getDbTfs <- function(scenicOptions)
#' {
#'   motifAnnotations <- getDbAnnotations(scenicOptions)
#'   
#'   allTFs <- sort(unique(motifAnnotations$TF))
#'   return(allTFs)
#' }
