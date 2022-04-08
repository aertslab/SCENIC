#' @title Object to store SCENIC settings 
#' @aliases loadInt show
#' @rdname ScenicOptions-class
#' @description
#' This class contains the options/settings for a run of SCENIC. 
#' Most SCENIC functions use this object as input instead of traditional arguments that need to be set individually. 
#' 
#' The object has three main slots:
#' \itemize{
#'    \item \code{@inputDatasetInfo}: Contains the information about the dataset to analyze: 
#'    dataset name ("datasetTitle", only for user information), 
#'    organism ("org", determines the motif databases to use), and 
#'    the files containing cell phenotype information ("cellInfo", "colVars", for plots. optional). 
#'    
#'    An overview of this slot can be obtained with \code{getDatasetInfo(scenicOptions)}.
#'    \item \code{@fileNames}: Contains the file names where the results are saved (\code{$output}: most relevant results, \code{$int}: intermediate files).
#'    
#'    Output file names can be obtained with \code{getOutName(scenicOptions)}. To load an intermediate file: \code{getIntName(scenicOptions)} and \code{regulons <- loadInt(scenicOptions, "aucell_regulons")}.
#'    \item \code{@settings}: Arguments for specific functions/steps:
#'    
#'    - General arguments ("verbose", "nCores"), and "seed" for AUCell rankings and t-SNEs.
#'    
#'    - \code{runSCENIC_1_coexNetwork2modules()}: "modules/weightThreshold" for the co-expression modules. 
#'    
#'    - \code{runSCENIC_2_createRegulons()}: RcisTarget databases ("dbs", "db" , "dbDir"). These are used in \code{runSCENIC_2_createRegulons()}, but the input expression matrix and GENIE3/GRNBoost regulators should be consistent.
#'    
#'    - \code{runSCENIC_3_scoreCells()}: "aucell/smallestPopPercent" for AUCell automatic thresholds.
#'    
#'    - t-SNEs: "defaultTsne/perpl", "defaultTsne/dims", "defaultTsne/aucType", "tSNE_filePrefix" (and "seed").
#'    
#'    The overview of this slot can be obtained with \code{getSettings(scenicOptions)}.
#' }
#'
#' In the current version there are not specific functions for setting a value. 
#' Follow the guidelines in the specific function if you need to modify a specific parameter.
#' 
#' @return
#' \itemize{
#' \item \code{initializeScenic()}: Creates the object. It also creates the folders where the results will be saved: 'int' and 'output'.
#' \item \code{show()}: Prints a summary of the object
#' \item \code{loadInt()}: Loads the selected "intermediate" file (normally from folder 'int/'). \code{getIntName(scenicOptions)} lists all possibilities (rownames: object name, fileName: file that will be loaded).
#' \item \code{getDatasetInfo()}, \code{getOutName()}: Shows the content of the corresponding slots.
#' }
#' @method show ScenicOptions
#' @method loadInt ScenicOptions
#' @examples
#' # Create object:
#' data(defaultDbNames)
#' scenicOptions <- initializeScenic(org="hgnc", datasetTitle="My dataset", 
#' dbDir="databases", dbs=defaultDbNames[["hgnc"]], nCores=4)
#' 
#' ### Accessor functions
#' # Get output file names:
#' getOutName(scenicOptions)  # Shows all
#' getOutName(scenicOptions, "s2_motifEnrichmentHtml") 
#' 
#' # Load intermediate files: 
#' getIntName(scenicOptions) # Shows all (use the rowname to subset)
#' regulons <- loadInt(scenicOptions, "aucell_regulons") # load the file
#' 
#' # other:
#' getDatasetInfo(scenicOptions)
#' getDatasetInfo(scenicOptions, "datasetTitle")
#' scenicOptions@inputDatasetInfo$datasetTitle <- "new title" # to assign new value
#' getSettings(scenicOptions, "defaultTsne/dims")  
#' 
#' tsneFileName(scenicOptions)
#' getDatabases(scenicOptions) 
#' dbVersion(getDatabases(scenicOptions))
#' @export ScenicOptions
ScenicOptions <- setClass(
  # Set the name for the class
  Class="ScenicOptions",
  
  slots = c(
    inputDatasetInfo="list",
    status="list",
    settings="list",
    fileNames="list"
  )
)

#' @importFrom R.utils capitalize
#' @rdname ScenicOptions-class
#' @aliases show,ScenicOptions-method
setMethod("show",
           signature="ScenicOptions",
           definition = function(object) {
            
             message <- paste("SCENIC settings object \n",
                              "\t Dataset title: ",getDatasetInfo(scenicOptions, "datasetTitle"), "\n",
                              "\t Organism: ", getDatasetInfo(scenicOptions, "org"), "\n", sep="")
             message <- paste(message,
                              "\t Cell info and colvars files: ",
                              getDatasetInfo(scenicOptions, "cellInfo"),", ",
                              getDatasetInfo(scenicOptions, "colVars"), "\n", sep="") 

            if(is.null(object@settings$dbs)){
              message <- paste(message,
                               "\nInitialise DB\n", sep="") 
            }else{
              message <- paste(message,
                               "\t Databases for this analysis: ",
                               paste0(getDatabases(scenicOptions), collapse=", "), "\n", sep="") 
            }
            
            message <- paste0(message,
                              "Current status: ",
                             "(",getStatus(scenicOptions, asID=TRUE),") " ,
                             getStatus(scenicOptions)) 
             
            cat(message) 
          }
)

##### Get dataset info ----
#' @name getDatasetInfo
#' @rdname ScenicOptions-class
#' @export getDatasetInfo 
setGeneric(name="getDatasetInfo",
           def=function(object, ...) standardGeneric("getDatasetInfo"))

#' @rdname ScenicOptions-class
#' @aliases getDatasetInfo,ScenicOptions-method
#' @exportMethod getDatasetInfo
setMethod("getDatasetInfo",
          signature="ScenicOptions",
          definition = function(object, slotName=NULL)
          {
            if(is.null(slotName))
            {
              print(do.call(rbind, object@inputDatasetInfo)[,1,drop=FALSE]) #TODO
            }else{
              if(slotName=="org") return(object@inputDatasetInfo$org)
              if(slotName=="datasetTitle") return(object@inputDatasetInfo$datasetTitle)
              if(slotName=="cellInfo") return(object@inputDatasetInfo$cellInfo[1])
              if(slotName=="colVars") return(object@inputDatasetInfo$colVars[1])
            }
            invisible(NULL)
          })

##### Get the path to the databases ----
# setGeneric
#' @name getDatabases
#' @rdname ScenicOptions-class
#' @export getDatabases 
setGeneric(name="getDatabases",
           def=function(object, ...) standardGeneric("getDatabases"))

#' @rdname ScenicOptions-class
#' @aliases getDatabases,ScenicOptions-method
#' @exportMethod getDatabases
setMethod("getDatabases",
          signature="ScenicOptions",
          definition =  function(object)
          {
            fullPath <- sapply(getSettings(object, "dbs"),function(x) file.path(getSettings(object, "dbDir"), x))
            return(fullPath)
          })

##### Get status ----
#' @name getStatus
#' @rdname ScenicOptions-class
#' @export getStatus 
setGeneric(name="getStatus", def=function(object, ...) standardGeneric("getStatus"))

#' @rdname ScenicOptions-class
#' @aliases getStatus,ScenicOptions-method
#' @exportMethod getStatus
setMethod("getStatus",
          signature="ScenicOptions",
          definition = function(object, asID=FALSE)
          {
            ret <- object@status$values[object@status$current]
            
            if(asID) ret <- object@status$current
            
            return(ret)
          })

##### Get settings ----
#' @name getSettings
#' @rdname ScenicOptions-class
#' @export getSettings 
setGeneric(name="getSettings",
           def=function(object, ...) standardGeneric("getSettings"))

#' @rdname ScenicOptions-class
#' @aliases getSettings,ScenicOptions-method
#' @exportMethod getSettings
setMethod("getSettings",
          signature="ScenicOptions",
          definition = function(object, slotName=NULL)
          {
            if(is.null(slotName))
            {
              print(object@settings)
            }else{
              if(slotName=="verbose") return(object@settings$verbose)
              if(slotName=="nCores") return(object@settings$nCores)
              if(slotName=="seed") return(object@settings$seed)
              if(slotName=="dbs") return(object@settings$dbs)
              if(slotName=="db") return(object@settings$dbs)
              if(slotName=="dbDir") return(path.expand(object@settings$dbDir))
              if(slotName=="modules/weightThreshold") return(object@settings$modules$weightThreshold)
              if(slotName=="aucell/smallestPopPercent") return(object@settings$aucell$smallestPopPercent)
              if(slotName=="defaultTsne/perpl") return(object@settings$defaultTsne$perpl)
              if(slotName=="defaultTsne/dims") return(object@settings$defaultTsne$dims)
              if(slotName=="defaultTsne/aucType") return(object@settings$defaultTsne$aucType)
              if(slotName=="tSNE_filePrefix") return(object@settings$tSNE_filePrefix)
              
              if(slotName=="devType") return(object@settings$devType)
            }
            invisible(NULL)
          })


##### Get the file name for an output ----
# setGeneric
#' @name getOutName
#' @rdname ScenicOptions-class
#' @export getOutName 
setGeneric(name="getOutName",
           def=function(object, ...) standardGeneric("getOutName"))

#' @rdname ScenicOptions-class
#' @aliases getOutName,ScenicOptions-method
#' @exportMethod getOutName
setMethod("getOutName",
          signature="ScenicOptions",
          definition = function(object, out_type=NULL)
          {
            if(is.null(out_type)) return(object@fileNames$output)
            object@fileNames$output[out_type,1]
          })

##### Get the file name for an intermediate step ----
# setGeneric 
#' @name getIntName
#' @rdname ScenicOptions-class
#' @export getIntName 
setGeneric(name="getIntName",
           def=function(object, ...) standardGeneric("getIntName"))

#' @rdname ScenicOptions-class
#' @aliases getIntName,ScenicOptions-method
#' @exportMethod getIntName
setMethod("getIntName",
          signature="ScenicOptions",
          definition = function(object, int_type=NULL)
            {
                ret <- ""
                
                if(is.null(int_type)) ret <- object@fileNames$int
                else if(int_type %in% rownames(object@fileNames$int))  ret <- object@fileNames$int[int_type,1]
                  
                return(ret)
              })


##### Load item from scenicOptions ----
# setGeneric
# @method test data.frame
#' @name loadFile
#' @rdname ScenicOptions-class
#' @export loadFile 
setGeneric(name="loadFile",
           def=function(object,...) standardGeneric("loadFile"))

#' @rdname ScenicOptions-class
#' @aliases loadFile,ScenicOptions-method
#' @exportMethod loadFile
setMethod("loadFile",
          signature="ScenicOptions",
          definition = function(object, fileName, verbose=FALSE, ifNotExists=c("error", "null"), ...)
          {
            retObj <- NULL
            
            if(!file.exists(as.character(fileName)) || is.null(fileName)) {
              ifNotExists <- ifNotExists[1]
              if(ifNotExists=="error") stop("File '", fileName, "' does not exist.")
              if(ifNotExists=="null")
              {
                if(verbose) message(fileName, " does not exist.")
                return(retObj)
              }
            }else{
              if(verbose) message("Loading ", fileName)
              tryCatch({
                if(grepl(".txt$", fileName))
                {
                  retObj <- read.table(fileName, sep="\t", ...)
                }else if(grepl(".rdata$", tolower(fileName))) {
                  retObj <- eval(as.name(load(fileName, verbose=TRUE)))
                }else {
                  retObj <- readRDS(fileName)
                }
              }, error=function(e){
                e$message=paste(e$message, "[file: ",fileName,"]")
                stop(e)
              })
              return(retObj)
            }
          }
)

##### Load intermediate result ----
# setGeneric
# @method test data.frame
#' @name loadInt
#' @rdname ScenicOptions-class
#' @export loadInt 
setGeneric(name="loadInt",
           def=function(object,...) standardGeneric("loadInt"))

#' @rdname ScenicOptions-class
#' @aliases loadInt,ScenicOptions-method
#' @exportMethod loadInt
setMethod("loadInt",
          signature="ScenicOptions",
          definition = function(object, int_type=NULL, ...)
          {
            if(is.null(int_type)) getIntName(object)
            else {
                loadFile(object, fileName=getIntName(object, int_type) , ...)
            }
          }
)


################## Functions ######################################
#' @rdname ScenicOptions-class
#' @export 
initializeScenic <- function(org=NULL, dbDir="databases", dbs=NULL, datasetTitle="", nCores=4, dbIndexCol='features')
{
  inputDataset<- list(
    org=org,
    datasetTitle=datasetTitle,
    cellInfo=c("int/cellInfo.Rds", NA),
    colVars=c("int/colVars.Rds", NA),
    int_01=c("int/cellColorNgenes.Rds", NA)
  )
  
  ## Load databases
  if(!org %in% c("mgi","hgnc", "dmel")) stop("'org' should be one of: mgi, hgnc, dmel.") 
  defaultDBs <- FALSE
  if(is.null(dbs))
  {
    data(defaultDbNames)
    dbs <- defaultDbNames[[org]]
  }
  dbsFound <- unlist(unname(lapply(dbs,function(x) setNames(file.exists(file.path(dbDir, x)), unname(file.path(dbDir, x))))))
  if(any(!dbsFound)) 
  {
    stop("The following RcisTarget databases were not found: ",
            paste(paste0("\n- ", names(dbsFound[which(!dbsFound)])), collapse=" "), 
         "\nMake sure the arguments 'dbDir' and 'dbs' are correct.")
    dbs <- NULL
  }else
  {
    message("Motif databases selected: ", 
            paste(paste0("\n  ", dbs, collapse=" ")))  
  }
  loadAttempt <- sapply(dbs,function(x) dbLoadingAttempt(file.path(dbDir, x), indexCol=dbIndexCol))
  if(any(!loadAttempt)) warning("It was not possible to load the following databses; check whether they are downloaded correctly: \n",
                                paste(dbs[which(!loadAttempt)], collapse="\n"))
  
  db_mcVersion <- dbVersion(dbs)
  
  scenicSettings=list(
    dbs=dbs,
    dbDir=dbDir,
    db_mcVersion=db_mcVersion,
    db_annotFiles=NULL,
    verbose=TRUE,
    nCores=nCores,
    seed=123,
    devType="pdf",
    modules=list(weightThreshold=0.001),
    regulons=list(),
    aucell=list(smallestPopPercent=0.25),
    defaultTsne=list(dims=50,
                perpl=30,
                aucType="AUC"),
    tSNE_filePrefix="int/tSNE"
  )
  
  # Files (for output and intermediate files)
  scenicFiles <- list(
    output=c(
      s2_motifEnrichment="output/Step2_MotifEnrichment.tsv", 
      s2_motifEnrichmentHtml="output/Step2_MotifEnrichment_preview.html", 
      s2_regulonTargetsInfo="output/Step2_regulonTargetsInfo.tsv", 
      s3_AUCheatmap="output/Step3_RegulonActivity_heatmap",
      s3_AUCtSNE_colAct="output/Step3_RegulonActivity_tSNE_colByActivity",
      s3_AUCtSNE_colProps="output/Step3_RegulonActivity_tSNE_colByCellProps",
      s4_boxplotBinaryActivity="output/Step4_BoxplotActiveCellsRegulon",
      s4_binaryActivityHeatmap="output/Step4_BinaryRegulonActivity_Heatmap_",
      s4_binarytSNE_colAct="output/Step4_BinaryRegulonActivity_tSNE_colByActivity",
      s4_binarytSNE_colProps="output/Step4_BinaryRegulonActivity_tSNE_colByCellProps",
      loomFile="output/scenic.loom"
    ),
    
    int=list(
      genesKept=c("int/1.1_genesKept.Rds", TRUE),
      corrMat=c("int/1.2_corrMat.Rds", TRUE),
      genie3wm=c("int/1.3_GENIE3_weightMatrix.Rds", FALSE),
      genie3ll=c("int/1.4_GENIE3_linkList.Rds", TRUE),
      genie3weighPlot=c("int/1.5_weightPlot", TRUE),   #TODO Dev
      tfModules_asDF=c("int/1.6_tfModules_asDF.Rds", TRUE), 
      
      tfModules_forEnrichment=c("int/2.1_tfModules_forMotifEnrichmet.Rds", FALSE),
      motifs_AUC=c("int/2.2_motifs_AUC.Rds", FALSE),
      motifEnrichment_full=c("int/2.3_motifEnrichment.Rds", FALSE),
      motifEnrichment_selfMotifs_wGenes=c("int/2.4_motifEnrichment_selfMotifs_wGenes.Rds", FALSE),
      regulonTargetsInfo=c("int/2.5_regulonTargetsInfo.Rds", NA),
      regulons=c("int/2.6_regulons_asGeneSet.Rds", NA),
      regulons_incidMat=c("int/2.6_regulons_asIncidMat.Rds", NA), # TODO check Keep?
      
      aucell_regulons=c("int/3.1_regulons_forAUCell.Rds", NA),
      aucell_genesStatsPlot=c("int/3.2_aucellGenesStats", NA),
      aucell_rankings=c("int/3.3_aucellRankings.Rds", NA),
      aucell_regulonAUC=c("int/3.4_regulonAUC.Rds", NA),
      # aucell_tsneAucPrefix=c("int/3.4_tsneRegulonActivity", NA),
      aucell_thresholds=c("int/3.5_AUCellThresholds.Rds", NA),
      aucell_thresholdsTxt=c("int/3.5_AUCellThresholds_Info.tsv", NA),
      
      aucell_binary_full=c("int/4.1_binaryRegulonActivity.Rds", NA),
      aucell_binary_nonDupl=c("int/4.2_binaryRegulonActivity_nonDupl.Rds", NA),
      aucell_regulonSelection=c("int/4.3_regulonSelections.Rds", NA),
      aucell_binaryRegulonOrder=c("int/4.4_binaryRegulonOrder.Rds", NA)
      # aucell_tsneBinaryPrefix=c("int/4.5_tsneRegulonActivity", NA),
    )
  )
  scenicFiles$output <- cbind(fileName=scenicFiles$output)
  scenicFiles$int <- do.call(rbind,scenicFiles$int) #as.data.frame()
  colnames(scenicFiles$int) <- c("fileName", "keep")
  scenicFiles$int <- scenicFiles$int[,"fileName", drop=FALSE]  # TODO: dedice if used
  
  dir.create("int", showWarnings=FALSE)
  dir.create("output", showWarnings=FALSE)
  object <- new("ScenicOptions",
                inputDatasetInfo=inputDataset,
                status=list(current=0, values=c("0"="SCENIC initialized", 
                                                "1"="Co-expression modules", 
                                                "2"="Regulons", 
                                                "3"="Cells scored", 
                                                "4"="SCENIC run completed")),
                settings=scenicSettings,
                fileNames=scenicFiles)
  
  
  ## Check if motif annotation and rankings potentially match
  motifAnnot <- getDbAnnotations(object)
  featuresWithAnnot <- checkAnnots(object, motifAnnot)
  if(any(featuresWithAnnot == 0)) message("Missing annotations for: \n", paste("\t", names(which(featuresWithAnnot==0))))
  
  ## Return
  return(object)
}

#' @rdname ScenicOptions-class
#' @export 
dbVersion <- function(dbs)
{
  dbVersion <- NULL
  if(all(grepl(".mc9nr.", dbs))) dbVersion <- "v9"
  if(all(grepl(".mc8nr.", dbs))) dbVersion <- "v8"
  return(dbVersion)
}

#' @rdname ScenicOptions-class
#' @export 
dbLoadingAttempt <- function(dbFilePath, indexCol='features'){
    ret <- FALSE
    ret <- tryCatch({
    rf <- arrow::ReadableFile$create(dbFilePath)
    fr <- arrow::FeatherReader$create(rf)
    fr$version
    genesInDb <- names(fr)
    randomCol <- sample(genesInDb, 1)
    fr$Read(randomCol)
    rnk <- RcisTarget::importRankings(dbFilePath, columns=randomCol, indexCol=indexCol)
    return(TRUE)
  }
  , error=function(e){
    print(e$message)
    return(FALSE)
  }
  )

  return(ret)
}

#' @rdname ScenicOptions-class
#' @export 
checkAnnots <- function(object, motifAnnot)
{
  allFeaturesInAnnot <- unlist(motifAnnot[,1]) # motif or track
  featuresWithAnnot <-  lapply(getDatabases(object), function(dbFile) 
  {
    rnktype = "features"	#TODO: add as option for custom dbs
    nRnks <- getRanking(RcisTarget::importRankings(dbFile, columns = rnktype))
    nRnks <- dplyr::pull(nRnks, rnktype)

    length(intersect(allFeaturesInAnnot,nRnks))/length(unique(c(allFeaturesInAnnot,nRnks)))
  })
  return(featuresWithAnnot)
}

