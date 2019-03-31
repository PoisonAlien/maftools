#' Plots decomposed mutational signatures
#'
#' @description Takes results from \code{\link{extractSignatures}} and plots decomposed mutational signatures as a barplot.
#'
#' @param nmfRes results from \code{\link{extractSignatures}}
#' @param contributions If TRUE plots contribution of signatures in each sample.
#' @param color colors for each Ti/Tv conversion class. Default NULL
#' @param patient_order User defined ordering of samples. Default NULL.
#' @param title_size size of title. Default 1.3
#' @param axis_lwd axis width. Default 2.
#' @param font_size font size. Default 1.2
#' @param show_title Default TRUE
#' @param show_barcodes Default FALSE
#' @param yaxisLim Default 0.3. If NA autoscales.
#' @param ... further plot options passed to \code{\link{barplot}}
#' @return Nothing
#' @seealso \code{\link{trinucleotideMatrix}} \code{\link{plotSignatures}}
#' @export
#'
plotSignatures = function(nmfRes = NULL, contributions = FALSE, color = NULL, patient_order = NULL,
                          font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3, ...){

  conv.mat.nmf.signatures = nmfRes$signatures
  contrib = nmfRes$contributions
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
      lo = layout(mat = matrix(data = c(1, 2), nrow = 2), heights = c(6, 2))
      par(mar = c(6, 4, 2, 1))
      b = barplot(contrib, axes = FALSE, horiz = FALSE, col = cols, border = NA, names.arg = rep("", ncol(contrib)))
      axis(side = 1, at = b, labels = colnames(contrib), lwd = 2, cex.axis = font_size,
           las = 2, line = 0.2, hadj = 0.8, font = 1, tick = FALSE)
      axis(side = 2, at = seq(0, 1, 0.25), lwd = 3, font = 1, las = 2, cex.axis = 0.9)
      mtext(text = "Signature exposures", side = 2, font = 1, cex = 1, line = 2.8)
      plot.new()
      par(mar = c(2, 3, 0, 0))
      legend(x = "left", legend = rownames(contrib), col = cols[1:nrow(contrib)],
             border = NA, bty = "n", pch = 15, xpd = TRUE, ncol = 1,
             cex = 1.2, pt.cex = 1.5, horiz = TRUE)
    }else{
      lo = layout(mat = matrix(data = c(1, 2), nrow = 2), heights = c(6, 2))
      par(mar = c(3, 4, 2, 1))
      b = barplot(contrib, axes = FALSE, horiz = FALSE, col = cols, border = NA, names.arg = rep("", ncol(contrib)))
      axis(side = 2, at = seq(0, 1, 0.25), lwd = 3, font = 1, las = 2, cex.axis = 0.9)
      mtext(text = "Signature exposure", side = 2, font = 1, cex = 1, line = 2.8)
      plot.new()
      par(mar = c(2, 3, 0, 0))
      legend(x = "left", legend = rownames(contrib), col = cols[1:nrow(contrib)],
             border = NA, bty = "n", pch = 15, xpd = TRUE, ncol = 1,
             cex = 1.2, pt.cex = 1.5, horiz = TRUE)
    }
  }else{
    aetiology = structure(list(aetiology = c("spontaneous deamination of 5-methylcytosine",
                                             "APOBEC Cytidine Deaminase (C>T)", "defects in DNA-DSB repair by HR",
                                             "exposure to tobacco (smoking) mutagens", "Unknown", "defective DNA mismatch repair",
                                             "UV exposure", "Unknown", "defects in polymerase-eta", "defects in polymerase POLE",
                                             "exposure to alkylating agents", "Unknown", "APOBEC Cytidine Deaminase (C>G)",
                                             "Unknown", "defective DNA mismatch repair", "Unknown", "Unknown",
                                             "Unknown", "Unknown", "defective DNA mismatch repair", "unknown",
                                             "exposure to aristolochic acid", "Unknown", "exposures to aflatoxin",
                                             "Unknown", "defective DNA mismatch repair", "Unknown", "Unknown",
                                             "exposure to tobacco (chewing) mutagens", "Unknown")), .Names = "aetiology", row.names = c("Signature_1",
                                                                                                                                        "Signature_2", "Signature_3", "Signature_4", "Signature_5", "Signature_6",
                                                                                                                                        "Signature_7", "Signature_8", "Signature_9", "Signature_10",
                                                                                                                                        "Signature_11", "Signature_12", "Signature_13", "Signature_14",
                                                                                                                                        "Signature_15", "Signature_16", "Signature_17", "Signature_18",
                                                                                                                                        "Signature_19", "Signature_20", "Signature_21", "Signature_22",
                                                                                                                                        "Signature_23", "Signature_24", "Signature_25", "Signature_26",
                                                                                                                                        "Signature_27", "Signature_28", "Signature_29", "Signature_30"
                                             ), class = "data.frame")

    plotData = as.data.frame(t(conv.mat.nmf.signatures))
    nsigs = nrow(plotData)

    if(is.null(color)){
      #color = c("blue","black","red","gray","green","maroon")
      color = c('coral4', 'lightcyan4', 'deeppink3', 'lightsalmon1', 'forestgreen', 'cornflowerblue')
    }
    colors = rep(color, each=16)

    par(mfrow = c(nsigs,1),oma = c(5,4,0,0) + 0.1, mar = c(0,0,2.5,0) + 0.1, las=1, tcl=-.25, font.main=4, xpd = NA)

    for(i in 1:nsigs){
      ae.sig = names(which(coSineMat[i,] == max(coSineMat[i,])))
      ae = as.character(aetiology[ae.sig,])
      #ae = paste0("Aetiology: ", ae, " \n cosine-similarity: ", max(coSineMat[i,]))
      ae = paste0(ae.sig, " like; cosine-similarity: ", round(max(coSineMat[i,]), digits = 3), " \n Aetiology: ", ae)
      d = as.matrix(plotData[i,])
      if(is.na(yaxisLim)){
        bh = ceiling(max(d, na.rm = TRUE) * 10)/10 #Bar height
      }else{
        bh = 0.3
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
