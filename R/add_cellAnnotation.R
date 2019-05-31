
#' @import SingleCellExperiment
#' @import SingleCellExperiment
load_as_sce <- function(loomPath)
{
  loom <- open_loom(loomPath, mode="r")
  dgem <- get_dgem(loom)
  cellAnnot <- get_cellAnnotation(loom)
  embeddings <- get_embeddings(loom)
  close_loom(loom)

  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=dgem),
                                                    colData=data.frame(cellAnnot[colnames(dgem),, drop=FALSE]),
                                                    reducedDims=S4Vectors::SimpleList(embeddings))

  return(sce)
}


add_cellAnnotation <- function(loom, cellAnnotation)
{
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}
