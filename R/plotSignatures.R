#' Plots decomposed mutational signatures.
#'
#' @description Plots decomposed mutational signatures as a barplot.
#'
#' @param nmfRes results from \code{\link{extractSignatures}}
#' @param contributions If TRUE plots contribution of signatures in each sample.
#' @param color colors for each Ti/Tv conversion class. Default NULL
#' @param ... further plot options passed to \code{\link{barplot}}
#' @return ggplot object if contributions is TRUE
#' @seealso \code{\link{trinucleotideMatrix}}
#' @export
plotSignatures = function(nmfRes = NULL, contributions = FALSE, color = NULL, ...){

  conv.mat.nmf.signatures = nmfRes$signatures
  contrib = nmfRes$contributions

  if(contributions){
    #     if(nrow(contrib) == 1){
    #       stop('Cannot plot contriubutions for single signature')
    #     }
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

    if(is.null(color)){
      #color = c("blue","black","red","gray","green","maroon")
      color = c('coral4', 'lightcyan4', 'deeppink3', 'lightsalmon1', 'forestgreen', 'cornflowerblue')
    }
    colors = rep(color, each=16)

    for(i in 1:nsigs){
      d = as.matrix(plotData[i,])
      barplot(d, xaxt = "n", yaxt = "n", border = FALSE, col = colors, beside = TRUE, ylim = c(-0.1, 0.2), ...)
      axis(side = 2,at = c(0,0.05,0.1,0.15,0.2),labels = c(0,0.05,0.1,0.15,0.2), pos = c(0,0.05,0.1,0.15,0.2), cex.axis=0.6, las = 2)
      abline(h = c(0.05,0.1,0.15,0.2,0.25),lty=2,lwd=0.3, col = 'gray70')
      rect(xleft = seq(0, 192, 32), ybottom = -0.05, xright = 192, ytop = -0.02, col = color, border = 'gray70')
      text(labels = c("C>A","C>G","C>T","T>A","T>C","T>G"),y = rep(-0.08,6),x = seq(0, 192, 32)[2:7]-16, cex = 0.6)
    }
  }
}
