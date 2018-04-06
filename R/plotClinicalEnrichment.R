#' Plots results from clinicalEnrichment analysis
#' @param enrich_res results from \code{\link{clinicalEnrichment}} or \code{\link{signatureEnrichment}}
#' @param pVal Default 0.05
#' @param cols named vector of colors for factor in a clinical feature. Default NULL
#' @param annoFontSize cex for annotation font size. Default 0.8
#' @param geneFontSize cex for gene font size. Default 0.8
#' @param legendFontSize cex for legend font size. Default 0.8
#' @param showTitle Default TRUE
#' @return returns nothing.
#' @export
#' @seealso \code{\link{clinicalEnrichment}} \code{\link{signatureEnrichment}}
#'
plotEnrichmentResults = function(enrich_res, pVal = 0.05, cols = NULL, annoFontSize = 0.8, geneFontSize = 0.8, legendFontSize = 0.8, showTitle = TRUE){

  res = enrich_res$groupwise_comparision

  plot.dat = data.table::data.table(
    Hugo_Symbol = res$Hugo_Symbol,
    g1_muts = as.numeric(sapply(strsplit(x = res$n_mutated_group1, split = " of "), "[[", 1)),
    g1_tot = as.numeric(sapply(strsplit(x = res$n_mutated_group1, split = " of "), "[[", 2)),
    g2_muts = as.numeric(sapply(strsplit(x = res$n_mutated_group2, split = " of "), "[[", 1)),
    g2_tot = as.numeric(sapply(strsplit(x = res$n_mutated_group2, split = " of "), "[[", 2)),
    P_value = res$p_value, Group1 = res$Group1, Group2 = "Res"
  )

  plot.dat$g1_muts_fract = apply(plot.dat, 1, function(x) round(as.numeric(x[2])/as.numeric(x[3]), digits = 3))
  plot.dat$g2_muts_fract = apply(plot.dat, 1, function(x) round(as.numeric(x[4])/as.numeric(x[5]), digits = 3))

  plot.dat[,g1_title := paste0(g1_muts, "/", g1_tot)]
  plot.dat[,g2_title := paste0(g2_muts, "/", g2_tot)]

  plot.dat = plot.dat[P_value < pVal]

  if(is.null(cols)){
    cols = RColorBrewer::brewer.pal(n = 9,name = 'Spectral')
    names(cols) = as.character(enrich_res$cf_sizes$cf)
  }else{
    names(cols) = as.character(enrich_res$cf_sizes$cf)
  }

  bar.cols = cols[plot.dat[,Group1]]
  legend.cols = cols[unique(plot.dat[,Group1])]

  # from: https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
  add_legend <- function(...) {
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
                mar=c(0, 0, 0, 0), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend(...)
  }

  par(bty="n", mgp = c(0.5,0.5,0), las=1, tcl=-.25, font.main=4, xpd=TRUE, mar = c(4,3,3.5,5)) #
  #plot(c(0, nrow(plot.dat)+15),c(0,0),xlim=c(0.5,33),las=2, ylim=c(-1,1),xlab="", xaxt="n", type="l")
  b = barplot(height = plot.dat$g1_muts_fract, ylim = c(-1, 1), axes = FALSE, border = 0.1, col = bar.cols)
  text(b, plot.dat$g1_muts_fract+0.03 , plot.dat$g1_title ,cex = annoFontSize, las = 2, srt=90, adj=0, xpd=TRUE, font = 2)
  axis(side = 2, at = seq(-1, 1, 0.25), labels = c(rev(seq(0, 1, 0.25)), seq(0, 1, 0.25)[2:5]),
       lwd = 3, font.axis = 2, cex = 1.5, font = 2)
  b2 = barplot(height = -plot.dat$g2_muts_fract, add = TRUE, axes = FALSE, border = 0.1)
  text(b, -plot.dat$g2_muts_fract-0.03 , plot.dat$g2_title ,cex = annoFontSize, las = 2, srt=90, adj=1, xpd=TRUE, font = 2)
  text(b, -1 , plot.dat$Hugo_Symbol ,cex = geneFontSize, las = 2, srt=45, adj=1, xpd=TRUE, font = 4)
  #par(xpd = T, mar = par()$mar + c(0,0,0,7))
  add_legend("topright", pt.lwd = 2,
         legend = c(names(legend.cols), "Rest"), fill = c(legend.cols, "gray70"),
         bty = "n", cex = legendFontSize, border=NA, xpd = TRUE, text.font = 4)
  if(showTitle){
    title(main = enrich_res$clinicalFeature, adj = 0, cex.main = 1, outer = FALSE)
  }
}
