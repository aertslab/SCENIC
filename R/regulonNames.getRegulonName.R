# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

# @import
#'
#' @title getRegulonName
#' @description Gets the regulon name for a given TF (only returns the 'extended' regulons if no directly-annotated regulon is available)
#' @param TFs Transcription factor name
#' @param allRegulonNames List of regulon names (e.g. rownames(AUC...))
#' @return Named character vector
#' @details ...
#' @seealso List of vignettes included in the package: \code{vignette(package="SCENIC")}
#' @example
#' reguNames <- c("Dlx1 (103g)", "Dlx1_extended (190g)", "Olig2_extended (29g)", "Sox9 (17g)")
#' getRegulonName("Dlx1", reguNames)
#' getRegulonName("Olig2", reguNames)
#' @export
getRegulonName <- function(TFs, allRegulonNames)
{
  ret <- setNames(sapply(TFs, function(TF) allRegulonNames[grep(paste(TF, " \\(", sep=""), allRegulonNames)]), TFs)
  ret <- c(ret, setNames(sapply(TFs, function(TF) allRegulonNames[grep(paste(TF, "_extended", sep=""), allRegulonNames)]), TFs))
  ret <- unlist(ret)

  ret <- unlist(onlyNonDirectExtended(ret))

  ret
}

