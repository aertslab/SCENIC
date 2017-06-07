# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

# @import
#'
#' @title onlyNonDirectExtended
#' @description TO DO
#' @param regulonNames ...
#' @return ...
#' @details ...
#' @seealso List of vignettes included in the package: \code{vignette(package="SCENIC")}
# @example
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
