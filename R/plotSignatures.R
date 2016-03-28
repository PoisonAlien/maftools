#' Plots decomposed mutational signatures.
#'
#' @description Plots decomposed mutational signatures as a barplot.
#'
#' @param nmfRes results from \code{\link{extractSignatures}}
#' @export
plotSignatures = function(nmfRes = NULL){

  conv.mat.nmf.signatures = nmfRes$signatures

  plotData = as.data.frame(t(conv.mat.nmf.signatures))
  nsigs = nrow(plotData)
  par(mfrow = c(nsigs,1),oma = c(5,4,0,0) + 0.1,mar = c(0,0,1,1) + 0.1)
  color = rep(c("blue","black","red","gray","green","maroon"),each=16)

  for(i in 1:nsigs){
    barplot2(as.matrix(plotData[i,]),col=color,beside = T,axisnames=F,ylim=c(-0.1,0.3),cex.axis = 0.6,axes=F)
    axis(side = 2,at = c(0,0.05,0.1,0.15,0.2,0.25),labels = c(0,0.05,0.1,0.15,0.2,0.25),pos = c(0,0.05,0.1,0.15,0.2,0.25),cex.axis=0.8)
    abline(h = c(0.05,0.1,0.15,0.2,0.25),lty=2,lwd=0.3)
    text(labels = c("C>A","C>G","C>T","T>A","T>C","T>G"),y = rep(-0.05,6),x = c(20,50,80,110,140,170))
  }
}
