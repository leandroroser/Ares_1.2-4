
 
aresPlot <- function(output_aresCalc, title = NULL) {
  
  if(is.null(title)) {
    title <- ""
  }
  
  nrind <- nrow(output_aresCalc)   # number of individuals
  xl <- c(1, nrind)
  yl <- c( min(output_aresCalc[,2]), max(output_aresCalc[,3]) )  # axis limits
  plot(1:nrind,output_aresCalc[,1], type = "l", col = "black", 
       ylab = "Allelic Richness", xlab = "Number of Individuals",
       xlim = xl, ylim = yl)
  title(main = title)
  lines(1:nrind,output_aresCalc[,2], type = "l", col = "red")
  lines(1:nrind,output_aresCalc[,3], type = "l", col = "red")
}