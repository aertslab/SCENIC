# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

# @import
#'
#' @title getTF
#' @description Returns the TF associated to a regulon name (e.g. removes the #genes sufix)
#' @param regulonName Character vector containing regulon names
#' @return Character vector containing the TF (including the "_extended" sufix if appropiate)
#' @details ...
#' @seealso List of vignettes included in the package: \code{vignette(package="SCENIC")}
#' @example
#' reguNames <- c("Dlx1 (103g)", "Dlx1_extended (190g)", "Olig2_extended (29g)", "Sox9 (17g)")
#' getTF(reguNames)
#' @export
getTF <- function(regulonName)
{
  gsub( "\\s\\(\\d+g)", "",regulonName)
}
