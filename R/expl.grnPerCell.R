# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

# @import
#' @title grnPerCell
#' @description TO DO
#' @param regulonsPerGroup ...
#' @param binaryMat ...
#' @param minRegulonPerc ...
#' @return ...
#' @details ...
#' @seealso List of vignettes included in the package: \code{vignette(package="SCENIC")}
# @example
#' @export
grnPerCell <- function(regulonsPerGroup, binaryMat, minRegulonPerc=0.3)
{
  regulonsPerGroup <- regulonsPerGroup[which(lengths(regulonsPerGroup)>0)]
  countCellRegulons <- matrix(ncol=ncol(binaryMat), nrow=length(regulonsPerGroup))
  colnames(countCellRegulons) <- colnames(binaryMat)
  rownames(countCellRegulons) <- names(regulonsPerGroup)
  for(cType in names(regulonsPerGroup))
  {
    regs <- regulonsPerGroup[[cType]]
    countCellRegulons[cType, colnames(binaryMat)] <-
      signif(apply(binaryMat[regs,, drop=F], 2, sum)/length(regs),1)
  }
  # table(apply(countCellRegulons, 2, sum))
  # table(apply(countCellRegulons, 2, function(x) sum(x==1)))
  # which(apply(countCellRegulons, 2, function(x) sum(x==1))==2)

  ### Assign cell to the highest GRN
  cellGroup <- setNames(rep("NA", ncol(countCellRegulons)), colnames(countCellRegulons))
  maxVals <- apply(countCellRegulons, 2, function(x) rownames(countCellRegulons)[which(x==max(c(minRegulonPerc, x)))]) # which.max: returns unique
  onlyOne <- unlist(maxVals[which(lengths(maxVals)==1)])
  cellGroup[names(onlyOne)] <- onlyOne

  moreThanOne <- reshape2::melt(maxVals[which(lengths(maxVals)>1)])
  cellGroup[unique(moreThanOne[,2])] <- "moreThanOne"

  ret <- list(cellGroup=cellGroup, moreThanOne=moreThanOne)
  return(ret)
}
