# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

# @import
#'
#' @title getTF
#' @description TO DO
#' @param regulonName ...
#' @return ...
#' @details ...
#' @seealso List of vignettes included in the package: \code{vignette(package="SCENIC")}
# @example
#' @export
getTF <- function(regulonName)
{
  gsub( "\\s\\(\\d+g)", "",regulonName)
}
