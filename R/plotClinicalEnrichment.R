#' Plots results from clinicalEnrichment analysis
#' @param enrich_res results from \code{\link{clinicalEnrichment}} or \code{\link{signatureEnrichment}}
#' @param pVal Default 0.05
#' @param ORthr Default 1. Odds ratio threshold. >1 indicates positive enrichment in the group of interest.
#' @param featureLvls Plot results from the selected levels. Default NULL, plots all.
#' @param cols named vector of colors for factor in a clinical feature. Default NULL
#' @param annoFontSize cex for annotation font size. Default 0.8
#' @param geneFontSize cex for gene font size. Default 0.8
#' @param legendFontSize cex for legend font size. Default 0.8
#' @param showTitle Default TRUE
#' @return returns nothing.
#' @export
#' @seealso \code{\link{clinicalEnrichment}} \code{\link{signatureEnrichment}}
#'
plotEnrichmentResults = function(enrich_res, pVal = 0.05, ORthr = 1, featureLvls = NULL, cols = NULL, annoFontSize = 0.8, geneFontSize = 0.8, legendFontSize = 0.8, showTitle = TRUE){

  res = enrich_res$groupwise_comparision

  plot.dat = data.table::data.table(
    Hugo_Symbol = res$Hugo_Symbol,
    g1_muts = as.numeric(sapply(strsplit(x = res$n_mutated_group1, split = " of "), "[[", 1)),
    g1_tot = as.numeric(sapply(strsplit(x = res$n_mutated_group1, split = " of "), "[[", 2)),
    g2_muts = as.numeric(sapply(strsplit(x = res$n_mutated_group2, split = " of "), "[[", 1)),
    g2_tot = as.numeric(sapply(strsplit(x = res$n_mutated_group2, split = " of "), "[[", 2)),
    P_value = res$p_value, Group1 = res$Group1, Group2 = "Res", OR = res$OR
  )

  plot.dat = plot.dat[P_value < pVal][OR > ORthr][order(Group1, -OR)]

  if(nrow(plot.dat) < 1){
    stop(paste0("No significant associations found at p-value < ", pVal, " and OR < ", ORthr))
  }

  if(!is.null(featureLvls)){
    plot.dat = plot.dat[Group1 %in% featureLvls]
    if(nrow(plot.dat) < 1){
      stop(paste0("No significant results found for ", paste(featureLvls, collapse = ", ")))
    }
  }

  conf_int_g1 = lapply(1:nrow(plot.dat), function(i){
    #estimate_binconf(X = plot.dat[i,g1_muts], n = plot.dat[i, g1_tot], alpha = pVal)
    as.data.frame(binconf(x = plot.dat[i,g1_muts], n = plot.dat[i, g1_tot], alpha = pVal))
  })
  conf_int_g1 = data.table::rbindlist(l = conf_int_g1)

  conf_int_g2 = lapply(1:nrow(plot.dat), function(i){
    #estimate_binconf(X = plot.dat[i,g2_muts], n = plot.dat[i, g2_tot], alpha = pVal)
    as.data.frame(binconf(x = plot.dat[i,g2_muts], n = plot.dat[i, g2_tot], alpha = pVal))
  })
  conf_int_g2 = data.table::rbindlist(l = conf_int_g2)

  plot.dat$g1_muts_fract = apply(plot.dat, 1, function(x) round(as.numeric(x[2])/as.numeric(x[3]), digits = 3))
  plot.dat$g2_muts_fract = apply(plot.dat, 1, function(x) round(as.numeric(x[4])/as.numeric(x[5]), digits = 3))

  plot.dat[,g1_title := paste0(g1_muts, "/", g1_tot)]
  plot.dat[,g2_title := paste0(g2_muts, "/", g2_tot)]

  if(is.null(cols)){
    cols = c(RColorBrewer::brewer.pal(n = 9,name = 'Set1'),
             RColorBrewer::brewer.pal(n = 8,name = 'Dark2'),
             RColorBrewer::brewer.pal(n = 8,name = 'Accent'))
    cols = cols[1:length(unique(plot.dat$Group1))]
    names(cols) = as.character(unique(plot.dat$Group1))
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

  yl_max = max(rbind(conf_int_g1, conf_int_g2), na.rm = TRUE)
  if(yl_max <= 1){
    yl_max = 1
  }
  data.table::setDF(x = conf_int_g1)
  data.table::setDF(x = conf_int_g2)

  #return(plot.dat)

  par(bty="n", mgp = c(0.5,0.5,0), las=1, tcl=-.25, font.main=4, xpd=TRUE, mar = c(4,3,3.5,1)) #
  #plot(c(0, nrow(plot.dat)+15),c(0,0),xlim=c(0.5,33),las=2, ylim=c(-1,1),xlab="", xaxt="n", type="l")
  b = barplot(height = plot.dat$g1_muts_fract, ylim = c(-1.25, yl_max), axes = FALSE, border = 0.1, col = bar.cols)

  #text(b, plot.dat$g1_muts_fract+0.03 , plot.dat$g1_title ,cex = annoFontSize, las = 2, srt=90, adj=0, xpd=TRUE, font = 1)
  axis(side = 2, at = seq(-1, 1, 0.25), labels = c(rev(seq(0, 1, 0.25)), seq(0, 1, 0.25)[2:5]),
       lwd = 1.2, font.axis = 2, cex = 1.5, font = 1)

  for(i in 1:nrow(conf_int_g1)){
    segments(x0 = b[i, 1], y0 = conf_int_g1[i,2],
             x1 = b[i, 1], y1 = conf_int_g1[i,3],
             lwd = 1.5)
  }
  text(b, conf_int_g1$Upper+0.03 , plot.dat$g1_title ,cex = annoFontSize, las = 2, srt=90, adj=0, xpd=TRUE, font = 1)

  b2 = barplot(height = -plot.dat$g2_muts_fract, add = TRUE, axes = FALSE, border = 0.1)
  for(i in 1:nrow(conf_int_g2)){
    segments(x0 = b2[i, 1], y0 = -conf_int_g2[i,2],
             x1 = b2[i, 1], y1 = -conf_int_g2[i,3],
             lwd = 1.5)
  }
  text(b, -conf_int_g2$Upper-0.03 , plot.dat$g2_title ,cex = annoFontSize, las = 2, srt=90, adj=1, xpd=TRUE, font = 1)

  text(b, -0.75 , plot.dat$Hugo_Symbol ,cex = geneFontSize, las = 2, srt = 90, adj = 1, xpd = TRUE, font = 3)

  # b = as.data.frame(b)
  # b$Group = plot.dat$Group1
  # b = split(b, b$Group)
  # for(idx in seq_along(b)){
  #   bg = b[[idx]]
  #   rect(xleft = bg$V1[1], ybottom = 1, xright = bg$V1[length(bg$V1)], ytop = 1.1, col = cols[idx])
  # }

  #par(xpd = T, mar = par()$mar + c(0,0,0,7))
  if(length(legend.cols) <= 4){
    n_col = 1
  }else{
    n_col = (length(legend.cols) %/% 4)+1
  }
  legend(x = 0, y = -1.1, pt.lwd = 2, ncol = n_col,
         legend = c(names(legend.cols), "Rest"), fill = c(legend.cols, "gray70"),
         bty = "n", cex = legendFontSize, border=NA, xpd = TRUE, text.font = 3)
  if(showTitle){
    title(main = enrich_res$clinicalFeature, adj = 0, cex.main = 1, outer = FALSE, font.main = 1)
  }
}
