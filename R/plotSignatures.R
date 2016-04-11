#' Plots decomposed mutational signatures.
#'
#' @description Plots decomposed mutational signatures as a barplot.
#'
#' @param nmfRes results from \code{\link{extractSignatures}}
#' @param contributions If TRUE plots contribution of signatures in each sample.
#' @param ... further plot options passed to \code{\link{barplot2}}
#' @return ggplot object if contributions is TRUE
#' @importFrom gplots barplot2
#' @export
plotSignatures = function(nmfRes = NULL, contributions = FALSE, ...){

  conv.mat.nmf.signatures = nmfRes$signatures
  contrib = nmfRes$contributions

  if(contributions){
    if(nrow(contrib) == 1){
      stop('Cannot plot contriubutions for single signature')
    }
    contribt = t(contrib)
    #calculate sd
    contribt = cbind(contribt, sd = apply(contribt, 1, sd))
    contribt = contribt[order(contribt[,ncol(contribt)]),] #order according to standard deviation
    contrib = t(contribt[,1:(ncol(contribt)-1)])
    contrib.melt = data.table::melt(contrib)

    contrib.gg = ggplot(data = contrib.melt, aes(x = Var2, y = value, fill = Var1))+geom_bar(stat = 'identity')+
      theme(legend.position = 'bottom', axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.title = element_blank())+ylab('Signature contribution')+xlab('Sample')

    return(contrib.gg)

  }else{
    plotData = as.data.frame(t(conv.mat.nmf.signatures))
    nsigs = nrow(plotData)
    par(mfrow = c(nsigs,1),oma = c(5,4,0,0) + 0.1,mar = c(0,0,1,1) + 0.1)
    color = rep(c("blue","black","red","gray","green","maroon"),each=16)

    for(i in 1:nsigs){
      gplots::barplot2(as.matrix(plotData[i,]),col=color,beside = TRUE, axisnames = FALSE, ylim=c(-0.1,0.3), cex.axis = 0.6, axes=FALSE, ...)
      axis(side = 2,at = c(0,0.05,0.1,0.15,0.2,0.25),labels = c(0,0.05,0.1,0.15,0.2,0.25), pos = c(0,0.05,0.1,0.15,0.2,0.25), cex.axis=0.8)
      abline(h = c(0.05,0.1,0.15,0.2,0.25),lty=2,lwd=0.3)
      text(labels = c("C>A","C>G","C>T","T>A","T>C","T>G"),y = rep(-0.05,6),x = c(20,50,80,110,140,170))
    }
  }
}
