# @import
#'
#' @title ...
#' @description  TO DO
#' To avoid calculating thresholds, set thresholds to FALSE
#' @param tSNE ...
#' @param exprMat ...
#' @param regulonAUC ...
#' @param thresholds ...
#' @param cex ...
#' @param alphaOn ...
#' @param alphaOff ...
#' @param borderColor ...
#' @param offColor ...
#' @param plots ...
#' @return ...
#' @details ...
#' @seealso List of vignettes included in the package: \code{vignette(package="SCENIC")}
# @example
#' @export
# thresholds=NULL; cex=1; alphaOn=1; alphaOff=0.2;  offColor="lightgray"
# borderColor=adjustcolor("darkgrey", alpha=.3); plots=c("histogram", "binaryAUC", "AUC", "expression")
plot_aucTsne <- function(tSNE, exprMat, regulonAUC=NULL, thresholds=NULL, cex=1,
                         alphaOn=1, alphaOff=0.2,
                         borderColor=adjustcolor("darkgrey", alpha=.3),
                         offColor="lightgray",
                         plots=c("histogram", "binaryAUC", "AUC", "expression"))
{
  library(BiocGenerics)
  library(AUCell)

  if(is.logical(thresholds) && thresholds == FALSE) thresholds <- NA
  if(!is.null(thresholds))
  {
    # if it is a list... probably return from AUCell. Let's try...
    if(is.list(thresholds[1])) {
      thresholds <- sapply(thresholds, function(x) unname(x$aucThr$selected))
    }

    if(!is.null(names(thresholds)))
    {
      regulons <- rownames(regulonAUC)[which(rownames(regulonAUC) %in% names(thresholds))]
      regulonAUC <- regulonAUC[regulons,]
    }
    if(is.null(names(thresholds)) || length(thresholds)==1)
    {
      thresholds <- setNames(rep(thresholds, nrow(regulonAUC)), rownames(regulonAUC))
    }
  }

  if(!is.null(regulonAUC)){
    selectedRegulons <- rownames(regulonAUC)
  }else{
    selectedRegulons <- rownames(exprMat)
    plots <- "expression"
    warning("no AUC provided... only plotting expression")
  }

  colorPal_Expr <- grDevices::colorRampPalette(c("orange", "darkorange", "brown", "brown"))
  colorPal_AUC <- grDevices::colorRampPalette(c("pink1", "indianred","red3", "darkred"))

  cells_trhAssignment <- list()
  for (regulon in selectedRegulons)
  {
    # Histogram
    if(is.null(thresholds) && any(c("histogram", "binaryAUC") %in% plots))
    {
      set.seed(123)
      cells_trhAssignment <- c(cells_trhAssignment,
          AUCell.exploreThresholds(regulonAUC[regulon,], assignCells=TRUE,
                                   plotHist=("histogram" %in% tolower(plots))))
      thisTrheshold <- cells_trhAssignment[[regulon]]$aucThr$selected
      thisAssignment <- cells_trhAssignment[[regulon]]$assignment
    }

    if(!is.null(thresholds) && !is.na(thresholds)) {
      if("histogram" %in% tolower(plots))
      {
        hist(getAuc(regulonAUC)[regulon,], main=regulon, sub="AUC", col="#00609060", border="#0060f0", breaks=100)
        abline(v=thresholds[regulon], lwd=3, lty=2, col="darkorange")
      }
      thisTrheshold <- thresholds[regulon]
      if(is.matrix(regulonAUC)){
        regulonAUCmatrix <- regulonAUC
      }else{
        regulonAUCmatrix <- getAuc(regulonAUC)
      }
      thisAssignment <- names(which(regulonAUCmatrix[regulon,] > thisTrheshold))
      cells_trhAssignment[[regulon]] <- list("threshold"=thisTrheshold, "assignment"=thisAssignment)
    }

    ######### Plot 1 #############
    # Cells assigned at current threshold
    if(any(grepl("binary", tolower(plots))))
    {
      pointBg <-  setNames(rep(adjustcolor("white" , alpha=alphaOff), nrow(tSNE)), rownames(tSNE))
      pointBg[thisAssignment] <- adjustcolor("midnightblue" , alpha=alphaOn)# "royalblue4"

      pointBorder <- setNames(rep(borderColor, nrow(tSNE)), rownames(tSNE))

      plot(tSNE, pch=21, col=pointBorder, bg=pointBg,
           axes=FALSE, ylab="", cex=cex,
           main=regulon,  xlab=paste("Cells with AUC > ", signif(thisTrheshold, 2), sep=""))
      box(which = "plot", col="grey")
    }

    ######### Plot 2 #############
    # Regulon AUC
    if("auc" %in% tolower(plots))
    {
      pointBorder <- setNames(rep(borderColor, nrow(tSNE)), rownames(tSNE))

      if (sum(getAuc(regulonAUC)[regulon,]) > 0)
      {
        # On
        cellColorMod <- setNames(adjustcolor(colorPal_AUC(5), alpha=alphaOn)[as.numeric(cut((getAuc(regulonAUC)[regulon,]), breaks=5, right=F,include.lowest=T))], colnames(regulonAUC))
        lowLim <- as.numeric(gsub("]", "", strsplit(levels(cut(getAuc(regulonAUC)[regulon,], breaks=100))[1], ",")[[1]][2]))
        # Off
        cellColorMod[which(getAuc(regulonAUC)[regulon,] <= max(0, lowLim))] <- adjustcolor(offColor, alpha=alphaOff)

        cellColorMod <- cellColorMod[rownames(tSNE)]

      } else
      {
        cellColorMod <- setNames(rep(adjustcolor(offColor, alpha=alphaOff), length(rownames(tSNE))), rownames(tSNE))
        pointBorder <- setNames(rep(borderColor, nrow(tSNE)), rownames(tSNE))
      }

      plot(tSNE, pch=21, col=pointBorder, bg=cellColorMod,
           axes=FALSE, ylab="", cex=cex,
           main=regulon, xlab="Regulon activity (AUC)")
      box(which = "plot", col="grey")
    }

    ######### Plot 3 #############
    # TF expression
    if("expression" %in% tolower(plots))
    {
      pointBorder <- setNames(rep(borderColor, nrow(tSNE)), rownames(tSNE))

      gene <- gsub( "\\s\\(\\d+g)", "",regulon)
      gene <- gsub( "_extended", "", gene)
      cellColorGene <- setNames(adjustcolor(colorPal_Expr(5)[as.numeric(cut(as.numeric(exprMat[gene,]), breaks=5, right=F, include.lowest=T))], alpha=alphaOn), colnames(exprMat))
      cellColorGene[which(exprMat[gene,]==0)] <- adjustcolor(offColor, alpha=alphaOff) # transparent
      cellColorGene <- cellColorGene[rownames(tSNE)]

      plot(tSNE, pch=21, col=pointBorder, bg=cellColorGene,
           axes=FALSE, xlab="", ylab="", cex=cex,
           main=paste(gene, "expression"))
      box(which = "plot", col="grey")
    }
  }
  invisible(cells_trhAssignment)
}

