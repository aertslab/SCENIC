#' @title RSS
#' @description Calculates the regulon specificity score
#' @param AUC
#' @param cellAnnotation
#' @param cellTypes
#' @return Matrix with the regulon specificity scores
#' @seealso 
#' The RSS was first used by Suo et al. in: 
#' Revealing the Critical Regulators of Cell Identity in the Mouse Cell Atlas. 
#' Cell Reports (2018). doi: 10.1016/j.celrep.2018.10.045
#' @examples 
#' # TODO
#' @export
calcRSS <- function(AUC, cellAnnotation, cellTypes=NULL)
{
  if(any(is.na(cellAnnotation))) stop("NAs in annotation")
  if(any(class(AUC)=="aucellResults")) AUC <- getAUC(AUC)
  normAUC <- AUC/rowSums(AUC)
  if(is.null(cellTypes)) cellTypes <- unique(cellAnnotation)
  # 
  ctapply <- lapply
  if(require('BiocParallel')) ctapply <- bplapply
  
  rss <- ctapply(cellTypes, function(thisType)
    sapply(rownames(normAUC), function(thisRegulon)
    {
      pRegulon <- normAUC[thisRegulon,]
      pCellType <- as.numeric(cellAnnotation==thisType)
      pCellType <- pCellType/sum(pCellType)
      .calcRSS.oneRegulon(pRegulon, pCellType)
    })
  )
  rss <- do.call(cbind, rss)
  colnames(rss) <- cellTypes
  return(rss)
}

#' @title Plot RSS
#' @description Plots an overview of the regulon specificity score
#' @param rss Output of calcRSS()
#' @param labelsToDiscard Cell labels to ignore (i.e. do not plot). IMPORTANT: All cells in the analysis should be included when calculating the RSS.
#' @param zThreshold 
#' @param cluster_columns 
#' @param order_rows 
#' @param trh 
#' @param varName 
#' @param col.low 
#' @param col.mid 
#' @param col.high 
#' @param setName Gene-set or cell type name to plot with plotRSS_oneSet()
#' @param n Number of top regulons to label in plotRSS_oneSet(). Default: 5
#' @return Matrix with the regulon specificity scores
#' @examples 
#' # TODO
#' @export
plotRSS <- function(rss, labelsToDiscard=NULL, zThreshold=1,
                    cluster_columns=FALSE, order_rows=TRUE, trh=0.01, varName="cellType",
                    col.low="grey90", col.mid="darkolivegreen3", col.high="darkgreen",
                    revCol=FALSE)
{
  varSize="RSS"
  varCol="Z"
  if(revCol) {
    varSize="Z"
    varCol="RSS"
  }
  
  rssNorm <- scale(rss) # scale the full matrix...
  rssNorm[rssNorm < zThreshold] <- 0
  rssNorm <- rssNorm[,which(!colnames(rssNorm) %in% labelsToDiscard)] # remove after calculating...
  
  ## to get topic order (easier...)
  tmp <- .plotRSS_heatmap(rssNorm, trh=trh, cluster_columns=cluster_columns, order_rows=order_rows)
  rowOrder <- rev(tmp@row_names_param$labels)
  
  ## Dotplot
  rss.df <- reshape2::melt(rss)
  head(rss.df)
  colnames(rss.df) <- c("Topic", varName, "RSS")
  rssNorm.df <- reshape2::melt(rssNorm)
  colnames(rssNorm.df) <- c("Topic", varName, "Z")
  rss.df <- merge(rss.df, rssNorm.df)
  
  rss.df <- rss.df[which(rss.df$Z >= 1.5),]
  rss.df <- rss.df[which(!rss.df[,varName] %in% labelsToDiscard),] # remove after calculating...
  # dim(rss.df)
  
  rss.df[,"Topic"] <- factor(rss.df[,"Topic"], levels=rowOrder)
<<<<<<< HEAD
  p <- dotHeatmap(rss.df, 
=======
  p <- dotheatmap(rss.df, 
>>>>>>> 611f9c2f9f3102623721ea5077ed74d3c709f6a7
             var.x=varName, var.y="Topic", 
             var.size=varSize, min.size=.5, max.size=5,
             var.col=varCol, col.low=col.low, col.mid=col.mid, col.high=col.high)
  
  invisible(list(plot=p, df=rss.df, rowOrder=rowOrder))
}

#' @aliases plotRSS
#' @export 
plotRSS_oneSet <- function(rss, setName, n=5)
{
  library(ggplot2)
  library(ggrepel)
  
  rssThisType <- sort(rss[,setName], decreasing=TRUE)
  thisRss <- data.frame(regulon=names(rssThisType), rank=seq_along(rssThisType), rss=rssThisType)
  thisRss$regulon[(n+1):nrow(thisRss)] <- NA
  
  ggplot(thisRss, aes(x=rank, y=rss)) + 
    geom_point(color = "blue", size = 1) + 
    ggtitle(setName) + 
    geom_label_repel(aes(label = regulon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50',
                     na.rm=TRUE) +
    theme_classic()
}



## Internal functions:
.H <- function(pVect){
  pVect <- pVect[pVect>0] # /sum(pVect) ??
  - sum(pVect * log2(pVect))
}

# Jensen-Shannon Divergence (JSD)
calcJSD <- function(pRegulon, pCellType)
{
  (.H((pRegulon+pCellType)/2)) - ((.H(pRegulon)+.H(pCellType))/2)
}

# Regulon specificity score (RSS)
.calcRSS.oneRegulon <- function(pRegulon, pCellType)
{
  jsd <- calcJSD(pRegulon, pCellType)
  1 - sqrt(jsd)
}

.plotRSS_heatmap <- plotRSS_heatmap <- function(rss, trh=NULL, row_names_gp=gpar(fontsize=5), order_rows=TRUE, cluster_rows=FALSE, name="RSS", ...)
{
  if(is.null(trh)) trh <- signif(quantile(rss, p=.97),2)
  
  library(ComplexHeatmap)
  rssSubset <- rss[rowSums(rss > trh)>0,]
  rssSubset <- rssSubset[,colSums(rssSubset > trh)>0]
  message("Showing regulons and cell types with any RSS > ", trh, " (dim: ", nrow(rssSubset), "x", ncol(rssSubset),")")
  
  if(order_rows)
  {
    maxVal <- apply(rssSubset, 1, which.max)
    rss_ordered <- rssSubset[0,]
    for(i in 1:ncol(rssSubset))
    {
      tmp <- rssSubset[which(maxVal==i),,drop=F]
      tmp <- tmp[order(tmp[,i], decreasing=FALSE),,drop=F]
      rss_ordered <- rbind(rss_ordered, tmp)
    }
    rssSubset <- rss_ordered
    cluster_rows=FALSE
  }
  
  Heatmap(rssSubset, name=name, row_names_gp=row_names_gp, cluster_rows=cluster_rows, ...)
} 

