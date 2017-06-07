# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

# @import
#'
#' @title getRegulonName
#' @description TO DO
#' @param TFs ...
#' @param allRegulonNames ...
#' @return ...
#' @details ...
#' @seealso List of vignettes included in the package: \code{vignette(package="SCENIC")}
# @example
#' @export
getRegulonName <- function(TFs, allRegulonNames)
{
  ret <- setNames(sapply(TFs, function(TF) allRegulonNames[grep(paste(TF, " \\(", sep=""), allRegulonNames)]), TFs)
  ret <- c(ret, setNames(sapply(TFs, function(TF) allRegulonNames[grep(paste(TF, "_extended", sep=""), allRegulonNames)]), TFs))
  ret <- unlist(ret)

  ret <- unlist(onlyNonDirectExtended(ret))

  ret
}

