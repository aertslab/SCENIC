#' @title dotHeatmap
#' @description Plots a dot-heatmap for enrichment results
#' @param enrichmentDf Input data.frame
#' @param var.x Variable (column) for the X axis
#' @param var.y Variable (column) for the Y axis
#' @param var.col Variable (column) that will determine the color of the dots
#' @param var.size Variable (column) that will determine the dot size
#' @param col.low Lower value color
#' @param col.mid Mid value color
#' @param col.high High value color
#' @param min.size Minimum dot size
#' @param max.size Maximum dot size
#' @param ... Other arguments to pass to ggplot's theme()
#' @return A ggplot object
#' @examples 
#' # TODO
#' @export
dotHeatmap <- function (enrichmentDf,
                        var.x="Topic", var.y="ID", 
                        var.col="FC", col.low="dodgerblue", col.mid="floralwhite", col.high="brown1", 
                        var.size="p.adjust", min.size=1, max.size=8,
                        ...)
{
  require(data.table)
  require(ggplot2)
  
  colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, col.high))
  p <- ggplot(data=enrichmentDf, mapping=aes_string(x=var.x, y=var.y)) + 
    geom_point(mapping=aes_string(size=var.size, color=var.col)) +
    scale_radius(range=c(min.size, max.size)) +
    scale_colour_gradientn(colors=colorPal(10)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y=element_blank(), 
          axis.text.x=element_text(angle=90, hjust=1),
          ...)
  return(p)
}

# temporary- TODO:delete
#' @export
dotheatmap <- dotHeatmap