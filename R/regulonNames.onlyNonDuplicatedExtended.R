# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

# @import
#'
#' @title onlyNonDuplicatedExtended
#' @description Returns the regulon names filtering-out the "extended" regulons if there is a regulon based on high-confidence annotations
#' @param regulonNames Character vector containing the regulon names (e.g. rownames(AUC_))
#' @return Character vector
#' @details ...
#' @seealso List of vignettes included in the package: \code{vignette(package="SCENIC")}
#' @examples
#' reguNames <- c("Dlx1 (103g)", "Dlx1_extended (190g)", "Olig2_extended (29g)", "Sox9 (17g)")
#' onlyNonDuplicatedExtended(reguNames)
#' @export
onlyNonDuplicatedExtended <- function(regulonNames)
{
  regulonNames <- unname(regulonNames)
  tfs <- getTF(regulonNames)
  tfs <- gsub("_extended", "", tfs)
  splitRegulons <- split(regulonNames, tfs)[unique(tfs)]

  ret <- sapply(splitRegulons, function(x) 
  {
    split(x, grepl("_extended", x))[[1]] # False (direct) will be first
  })
  
  return(ret)
}
