# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

# @import
#'
#' @title onlyNonDirectExtended
#' @description Returns the regulon names filtering-out the "extended" regulons if there is a regulon based on directly annotated motifs
#' @param regulonNames Character vector containing the regulon names (e.g. rownames(AUC_))
#' @return Character vector
#' @details ...
#' @seealso List of vignettes included in the package: \code{vignette(package="SCENIC")}
#' @examples
#' reguNames <- c("Dlx1 (103g)", "Dlx1_extended (190g)", "Olig2_extended (29g)", "Sox9 (17g)")
#' onlyNonDirectExtended(reguNames)
#' @export
onlyNonDirectExtended <- function(regulonNames)
{
  regulonNames <- unname(regulonNames)
  tfs <- gsub( "\\s\\(\\d+g)", "", regulonNames)
  tfs <- gsub( "_extended", "", tfs)
  splitRegulons <- split(regulonNames, tfs)

  direct <- sapply(splitRegulons[which(lengths(splitRegulons)>1)], function(x) x[which(!grepl("_extended", x))])

  any <- unlist(splitRegulons[which(lengths(splitRegulons)==1)])
  ret <- c(direct, any)
  ret[order(names(ret))]
}
