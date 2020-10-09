.openDev <- function(fileName, devType, ...)
{
  if(devType=="pdf")
    pdf(paste0(fileName, ".pdf"), ...)
  
  if(devType=="png")
    png(paste0(fileName, ".png", type="cairo"), ...)
  
  if(devType=="cairo_pfd") # similar to Cairo::CairoPDF?
    grDevices::cairo_pdf(paste0(fileName, ".pdf"), ...)
}

.openDevHeatmap <- function(fileName, devType)
{
  if(devType!="pdf") 
  {
    if(devType=="png") .openDev(fileName=fileName, devType=devType, width=1200,height=1200)
    if(devType!="png") .openDev(fileName=fileName, devType=devType)
    fileName <- NA
  }else{
    fileName <- paste0(fileName,".pdf")
  }
  return(fileName)
}

.closeDevHeatmap <- function(devType)
{
  if(devType!="pdf") 
  {
    dev.off()
  }
}