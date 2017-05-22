#' Plots decomposed mutational signatures or APOBEC enrichment plot.
#'
#' @description If input is results from \code{\link{extractSignatures}} plots decomposed mutational signatures as a barplot. If input is results from \code{\link{trinucleotideMatrix}}
#' plots APOBEC enrichment plot.
#'
#' @param nmfRes results from \code{\link{extractSignatures}} or \code{\link{trinucleotideMatrix}}
#' @param contributions If TRUE plots contribution of signatures in each sample.
#' @param color colors for each Ti/Tv conversion class. Default NULL
#' @param ... further plot options passed to \code{\link{barplot}}
#' @return ggplot object if contributions is TRUE
#' @seealso \code{\link{trinucleotideMatrix}}
#' @export
#'
plotSignatures = function(nmfRes = NULL, contributions = FALSE, color = NULL, ...){

  if(length(nmfRes) == 2){
    sub.tbl <- nmfRes$APOBEC_scores
    sub.tbl$APOBEC_Enriched = factor(sub.tbl$APOBEC_Enriched, levels = c('yes', 'no')) #Set levels
    #yp = boxplot.stats(x = sub.tbl[,n_mutations])$stats #yaxis points and limits
    yp = pretty(x = c(1: max(sub.tbl[,n_mutations], na.rm = TRUE)))
    yp[length(yp)] = max(sub.tbl[,n_mutations], na.rm = TRUE)
    yp[1] = 1

    if(nrow(sub.tbl[!is.na(APOBEC_Enriched), mean(fraction_APOBEC_mutations), APOBEC_Enriched][APOBEC_Enriched %in% 'yes']) == 0){
      stop('None of the samples are enriched for APOBEC. Nothing to plot.')
    }

    pieDat = sub.tbl[!is.na(APOBEC_Enriched), mean(fraction_APOBEC_mutations), APOBEC_Enriched]
    pieDat[,nonApobec := 1 - V1]
    colnames(pieDat)[2] = 'Apobec'
    pieDat = data.table::melt(pieDat, id.vars = 'APOBEC_Enriched', drop = FALSE)
    pieDat[,title := paste0(variable, ' [', round(value, digits = 3), ']')]
    pieDat$title = gsub(pattern = '^Apobec', replacement = 'tCw', x = pieDat$title)
    pieDat$title = gsub(pattern = '^nonApobec', replacement = 'non-tCw', x = pieDat$title)

    layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE), widths=c(2, 3))
    par(bty="n", mgp = c(0.5,0.5,0), las=1, tcl=-.25, font.main=4, xpd=NA)

    pieCol  = c("#084594", "#9ECAE1")

    boxplot(n_mutations ~ APOBEC_Enriched, data = sub.tbl,  xaxt="n", boxwex=0.6, outline = TRUE, lty=1,
            outwex=0, staplewex=0, frame.plot = FALSE, col = c('maroon', 'royalblue'), yaxt = 'n',
            ylim = c(min(yp), max(yp)),
            outcol="gray70", outcex = 0.8, outpch  = 16)
    title(main = 'Mutation load between APOBEC enriched \n and non-APOBEC enriched samples', cex.main=0.9)

    axis(side = 1, at = c(1, 2), labels = na.omit(sub.tbl[,.N,APOBEC_Enriched])[,paste0('N=', N)], las = 1, tick = FALSE)
    axis(side = 2, at = yp, lwd = 1.8, las = 1)

    pie(x = pieDat[APOBEC_Enriched %in% 'yes', value], col = pieCol,
        border="white", radius = 0.95, cex.main=0.6, labels =  pieDat[APOBEC_Enriched %in% 'yes', title], clockwise = TRUE)
    symbols(0,0,circles=.4, inches=FALSE, col="white", bg="white", lty=0, add=TRUE)
    title(main = 'Average tCw mutations in \n APOBEC enriched samples', cex.main=0.9)

    pie(x = pieDat[APOBEC_Enriched %in% 'no', value], col = pieCol,
        border="white",  radius = 0.95, cex.main=1.33, labels =  pieDat[APOBEC_Enriched %in% 'no', title], clockwise = TRUE)
    symbols(0,0,circles=.4, inches=FALSE, col="white", bg="white", lty=0, add=TRUE)
    title(main = 'Average tCw mutations in \n non-APOBEC enriched samples', cex.main = 0.9)

  }else{
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
        axis(side = 2,at = c(0,0.05,0.1,0.15,0.2),labels = c(0,0.05,0.1,0.15,0.2), pos = c(0,0.05,0.1,0.15,0.2), las = 2)
        abline(h = c(0.05,0.1,0.15,0.2,0.25),lty=2,lwd=0.3, col = 'gray70')
        rect(xleft = seq(0, 192, 32), ybottom = -0.05, xright = 192, ytop = -0.02, col = color, border = 'gray70')
        text(labels = c("C>A","C>G","C>T","T>A","T>C","T>G"),y = rep(-0.08,6),x = seq(0, 192, 32)[2:7]-16, cex = 0.8)
      }
    }
  }
}
