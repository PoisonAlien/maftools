#' Draw forest plot for differences betweeen cohorts.
#'
#' @details Plots results from \code{link{mafCompare}} as a forest plot with x-axis as log10 converted odds ratio and differentially mutated genes on y-axis.
#' @param mafCompareRes results from \code{\link{mafCompare}}
#' @param pVal p-value threshold. Default 0.05.
#' @param fdr fdr threshold. Default NULL. If provided uses adjusted pvalues (fdr).
#' @param show can be either \code{stat} or \code{pval}
#' @param color vector of colors for cohorts. Default NULL.
#' @param geneFontSize Font size for gene symbols. Default 12
#' @param file basename for output file. Plot will saved to an output pdf.
#' @param width width of plot to be generated
#' @param height height of plot to be generated
#' @export
#' @return ggplot object of the plot.
#' @seealso \code{\link{mafCompare}}
#' @examples
#' ##Primary and Relapse APL
#' primary.apl <- system.file("extdata", "APL_primary.maf.gz", package = "maftools")
#' relapse.apl <- system.file("extdata", "APL_relapse.maf.gz", package = "maftools")
#' ##Read mafs
#' primary.apl <- read.maf(maf = primary.apl)
#' relapse.apl <- read.maf(maf = relapse.apl)
#' ##Perform analysis and draw forest plot.
#' pt.vs.rt <- mafCompare(m1 = primary.apl, m2 = relapse.apl, m1Name = 'Primary',
#' m2Name = 'Relapse', minMut = 5)
#' forestPlot(mafCompareRes = pt.vs.rt, show = 'stat')

forestPlot = function(mafCompareRes, pVal = 0.05, fdr = NULL, show = NULL, color = NULL, geneFontSize = 12, file = NULL, width = 5, height = 6){

  res = mafCompareRes$results

  if(is.null(fdr)){
    m.sigs = res[pval < pVal]
  }else{
    m.sigs = res[adjPval < fdr]
  }

  m1Name = mafCompareRes$SampleSummary[1, Cohort]
  m2Name = mafCompareRes$SampleSummary[2, Cohort]

  m1.sampleSize = mafCompareRes$SampleSummary[1, SampleSize]
  m2.sampleSize = mafCompareRes$SampleSummary[2, SampleSize]

  if(nrow(m.sigs) < 1){
    stop('No differetially mutated genes found !')
  }

  m.sigs = fisherCorrection(fc = m.sigs)

  m.sigs$Hugo_Symbol = factor(x = m.sigs$Hugo_Symbol, levels = rev(m.sigs$Hugo_Symbol))
  m.sigs[,log10OR := log10(or)]
  m.sigs$label = paste('pval: ',round(m.sigs$pval, digits = 5), sep = '')
  m.sigs$flow = ifelse(test = m.sigs$log10OR < 0, yes = m2Name, no = m1Name)
  m.sigs$statRight = paste(m2Name,':' , m.sigs[,3,with =FALSE][[1]], sep = '')
  m.sigs$statLeft = paste(m1Name,':' , m.sigs[,2,with =FALSE][[1]], sep = '')


  if(!is.null(show)){
    if(show == 'pval'){
      m.sigs$label = paste('pval: ',round(m.sigs$pval, digits = 5), sep = '')
    }else if(show == 'stat'){
      m.sigs$label = apply(m.sigs[,.(statLeft, statRight)], 1, paste, collapse = ' ; ')
    }else{
      stop('show can only be pval or stat!')
    }
  }

  lim = max(abs(c(log10(m.sigs$ci.up), log10(m.sigs$ci.low))))+1

  gg.fp = ggplot(data = m.sigs, aes(x = Hugo_Symbol, y = log10OR, label = label, color = flow))+geom_point(size = 3)+
    geom_errorbar(aes(ymin = log10(ci.low), ymax = log10(ci.up)), size = 0.5, width = 0.20)+
    coord_flip()+cowplot::theme_cowplot(font_size = 9, line_size = 1)+cowplot::background_grid(major = 'x')+
    geom_hline(yintercept = 0, linetype = 'dotted')+xlab('Gene')+ylab('log10 (Odds Ratio)')+
    theme(axis.line.y = element_blank(), axis.text.y = element_text(face = "bold", size = geneFontSize), axis.text.x = element_text(face = "bold", size = 12),
          axis.title.x = element_text(face = "bold", size = 12), axis.title.y = element_blank(),
          legend.position = 'bottom', legend.title = element_blank(), legend.text = element_text(face = "bold", size = 12),
          plot.title = element_text(face = "bold", size = 14))+
    ylim(-lim, lim)

  if(!is.null(color)){
    color = color
    names(color) = c(m1Name, m2Name)
    gg.fp = gg.fp+scale_colour_manual(values = color)
  }

  if(!is.null(show)){
    gg.fp = gg.fp+ggrepel::geom_label_repel(size = 2.5, nudge_x = 0.2, force = 10, show.legend = FALSE, label.size = 0.2)
  }

  title = paste(m2Name, ' (n = ', m2.sampleSize, ')', ' v/s ' , m1Name, ' (n = ' ,m1.sampleSize, ')', sep='')
  #gg.fp = gg.fp+annotate(geom = 'text', y = c(-2, 2), x = c(nrow(m.sigs)+2, nrow(m.sigs)+2), label = c(m1Name, m2Name))
  gg.fp = gg.fp+ggtitle(label = title)
  print(gg.fp)

  if(!is.null(file)){
    cowplot::save_plot(filename = paste(file, 'pdf', sep='.'),
                       plot = gg.fp, base_height = height, base_width = width,
                       paper = "special", bg  = "white")
  }

  return(gg.fp)
}
