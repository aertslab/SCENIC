
# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

# @import
#' @title selectCells
#' @description TO DO
#' @param tSNE_binary ...
#' @param binaryRegulonActivity ...
#' @param verbose ...
#' @return ...
#' @details ...
#' @seealso List of vignettes included in the package: \code{vignette(package="SCENIC")}
# @example
#' @export
selectCells <- function(tSNE_binary, binaryRegulonActivity=NULL, verbose=TRUE)
{
    selectedCells <- gatepoints::fhs(tSNE_binary, mark = TRUE)

    ret <- list()
    ret[["cells"]] <- as.character(selectedCells)
    if(!is.null(binaryRegulonActivity))
    {
        nCellsOn <- sort(rowSums(binaryRegulonActivity[,selectedCells]), decreasing=TRUE)
        nCellsOn <- nCellsOn[nCellsOn>0]
        nCellsOther <- rowSums(binaryRegulonActivity[,which(!colnames(binaryRegulonActivity) %in% selectedCells), drop=F])
        ret[["regulons"]] <- cbind(nCellsOn=nCellsOn,
                                   percentage=nCellsOn/length(selectedCells),
                                   nCellsOnOtherGroups=nCellsOther[names(nCellsOn)],
                                   diff=(nCellsOn-nCellsOther[names(nCellsOn)])/nCellsOn)
    }

    if(verbose){
        if(!is.null(binaryRegulonActivity))
        {
            message(length(selectedCells), " cells selected.\n",
                    sum(ret[["regulons"]][,"percentage"]>0.90)," regulons active in more than 90% of them:")
            print(ret[["regulons"]][which(ret[["regulons"]][,"percentage"]>0.90),, drop=FALSE])
        }else{
            message(length(selectedCells), " cells selected.\n")
        }
    }
    return(ret)
}
