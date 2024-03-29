#' Plots decomposed mutational signatures
#'
#' @description Takes results from \code{\link{extractSignatures}} and plots decomposed mutational signatures as a barplot.
#'
#' @param nmfRes results from \code{\link{extractSignatures}}
#' @param contributions If TRUE plots contribution of signatures in each sample.
#' @param absolute Whether to plot absolute contributions. Default FALSE.
#' @param color colors for each Ti/Tv conversion class. Default NULL
#' @param patient_order User defined ordering of samples. Default NULL.
#' @param title_size size of title. Default 1.3
#' @param axis_lwd axis width. Default 1.
#' @param font_size font size. Default 0.6
#' @param show_title If TRUE compares signatures to COSMIC signatures and prints them as title
#' @param sig_db Only applicable if show_title is TRUE. can be \code{legacy}, \code{SBS}, \code{SBS_v34}. Default \code{SBS_v34}
#' @param show_barcodes Default FALSE
#' @param yaxisLim Default NA.
#' @param ... further plot options passed to \code{\link{barplot}}
#' @return Nothing
#' @seealso \code{\link{trinucleotideMatrix}} \code{\link{plotSignatures}}
#' @export
#'
plotSignatures = function(nmfRes = NULL, contributions = FALSE, absolute = FALSE, color = NULL, patient_order = NULL,
                          font_size = 0.6, show_title = TRUE, sig_db = "SBS_v34", axis_lwd = 1, title_size = 0.9, show_barcodes = FALSE, yaxisLim = NA, ...){

  conv.mat.nmf.signatures = nmfRes$signatures
  if(absolute){
    contrib = nmfRes$contributions_abs
  }else{
    contrib = nmfRes$contributions
  }

  coSineMat = nmfRes$coSineSimMat

  if(contributions){
    contribt = t(contrib)
    #calculate sd
    if(!is.null(patient_order)){
      contribt = contribt[patient_order,] #order on user-specified ordering of the genomes
    }else{
      contribt = contribt[order(contribt[,ncol(contribt)]),] #order according to standard deviation
    }

    #contrib = t(contribt[,1:(ncol(contribt)-1)])
    contrib = t(contribt[,1:(ncol(contribt))])

    cols = RColorBrewer::brewer.pal(n = 8, name = 'Set2')

    if(show_barcodes){
      lo = graphics::layout(mat = matrix(data = c(1, 2), nrow = 2), heights = c(6, 2))
      par(mar = c(6, 4, 2, 1))
      b = barplot(contrib, axes = FALSE, horiz = FALSE, col = cols, border = NA, names.arg = rep("", ncol(contrib)))
      axis(side = 1, at = b, labels = colnames(contrib), lwd = 2, cex.axis = font_size,
           las = 2, line = 0.2, hadj = 0.8, font = 1, tick = FALSE)
      axis(side = 2, at = pretty(contrib), lwd = 1, font = 1, las = 2, cex.axis = 0.9)
      mtext(text = "Signature exposures", side = 2, font = 1, cex = 1, line = 2.8)
      plot.new()
      par(mar = c(2, 3, 0, 0))
      legend(x = "left", legend = rownames(contrib), col = cols[1:nrow(contrib)],
             border = NA, bty = "n", pch = 15, xpd = TRUE, ncol = 1,
             cex = 1.2, pt.cex = 1.5, horiz = TRUE)
    }else{
      lo = graphics::layout(mat = matrix(data = c(1, 2), nrow = 2), heights = c(6, 2))
      par(mar = c(3, 4, 2, 1))
      b = barplot(contrib, axes = FALSE, horiz = FALSE, col = cols, border = NA, names.arg = rep("", ncol(contrib)))
      axis(side = 2, at = pretty(contrib), lwd = 1, font = 1, las = 2, cex.axis = 0.9)
      mtext(text = "Signature exposure", side = 2, font = 1, cex = 1, line = 2.8)
      plot.new()
      par(mar = c(2, 3, 0, 0))
      legend(x = "left", legend = rownames(contrib), col = cols[1:nrow(contrib)],
             border = NA, bty = "n", pch = 15, xpd = TRUE, ncol = 1,
             cex = 1.2, pt.cex = 1.5, horiz = TRUE)
    }
  }else{
    comp_res = compareSignatures(nmfRes = nmfRes, verbose = FALSE, sig_db = sig_db)

    plotData = as.data.frame(t(conv.mat.nmf.signatures))
    nsigs = nrow(plotData)

    if(is.null(color)){
      #color = c("blue","black","red","gray","green","maroon")
      color = c('coral4', 'lightcyan4', 'deeppink3', 'lightsalmon1', 'forestgreen', 'cornflowerblue')
    }
    colors = rep(color, each=16)

    par(mfrow = c(nsigs,1),oma = c(5,4,0,0) + 0.1, mar = c(0,0,2.5,0) + 0.1, las=1, tcl=-.25, font.main=4, xpd = NA)

    for(i in 1:nsigs){
      ae.bm = comp_res$best_match[[i]][2]
      ae = comp_res$best_match[[i]][1]
      #ae = paste0("Aetiology: ", ae, " \n cosine-similarity: ", max(coSineMat[i,]))
      ae = paste0(ae.bm, " \n Aetiology: ", ae)
      d = as.matrix(plotData[i,])
      if(is.na(yaxisLim)){
        bh = ceiling(max(d, na.rm = TRUE) * 10)/10 #Bar height
      }else{
        bh = yaxisLim
      }

      barplot(d, xaxt = "n", yaxt = "n", col = colors, beside = TRUE, ylim = c(-0.1, bh),
              cex.main = 1, border = NA, font.axis = 2, font.lab = 2,
              adj = 0.25, ...)
      if(show_title){
        title(main = ae, cex.main = title_size, line = 0, font.main = 3)
      }

      #mtext(text = ae, side = 1, line = 2, font = 1, cex = 0.5, at = 0.3)
      axis(side = 2, at = seq(0, bh, 0.1),
           pos = -2, las = 2, lwd = axis_lwd, hadj = 1.1,
           font = 1, cex.axis = font_size)
      #abline(h = seq(0, 0.3, 0.1),lty=2,lwd=0.3, col = 'gray70')
      rect(xleft = seq(0, 192, 32), ybottom = -0.05, xright = 192, ytop = -0.02, col = color, border = 'gray70')
      if(i == nsigs){
        text(labels = c("C>A","C>G","C>T","T>A","T>C","T>G"),
             y = rep(-0.1,6),x = seq(0, 192, 32)[2:7]-16, cex = font_size,
             font = 1, font.lab = 2, pos = 1.2)
      }
    }
  }
}
