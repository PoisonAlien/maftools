#' Draw forest plot for differences betweeen cohorts.
#'
#' @details Plots results from \code{link{mafCompare}} as a forest plot with x-axis as log10 converted odds ratio and differentially mutated genes on y-axis.
#' @param mafCompareRes results from \code{\link{mafCompare}}
#' @param pVal p-value threshold. Default 0.05.
#' @param fdr fdr threshold. Default NULL. If provided uses adjusted pvalues (fdr).
#' @param color vector of colors for cohorts. Default NULL.
#' @param geneFontSize Font size for gene symbols. Default 1.2
#' @param titleSize font size for titles. Default 1.2
#' @param lineWidth line width for CI bars. Default 2.2
#' @export
#' @return Nothing
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
#' forestPlot(mafCompareRes = pt.vs.rt)

forestPlot = function(mafCompareRes, pVal = 0.05, fdr = NULL, color = NULL,
                      geneFontSize = 1.2, titleSize = 1.2, lineWidth = 2.2){

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
  m.sigs[,log10OR_high := log10(ci.up)]
  m.sigs[,log10OR_low := log10(ci.low)]
  #m.sigs$label = paste('pval: ',round(m.sigs$pval, digits = 5), sep = '')
  m.sigs$flow = ifelse(test = m.sigs$log10OR < 0, yes = m2Name, no = m1Name)
  m.sigs$statRight = paste(m2Name,':' , m.sigs[,3,with =FALSE][[1]], sep = '')
  m.sigs$statLeft = paste(m1Name,':' , m.sigs[,2,with =FALSE][[1]], sep = '')
  m.sigs = m.sigs[order(pval, decreasing = TRUE)]


  lim = max(abs(c(m.sigs$log10OR_high[!is.infinite(m.sigs$log10OR_high)],
                  m.sigs$log10OR_low[!is.infinite(m.sigs$log10OR_low)])))+1
  lim = round(lim, digits = 2)

  m.sigs$log10OR = ifelse(test = is.infinite(m.sigs$log10OR),
                       yes = ifelse(test = m.sigs$log10OR > 0, yes = lim, no = -lim),
                       no = m.sigs$log10OR)

  m.sigs$log10OR_high = ifelse(test = is.infinite(m.sigs$log10OR_high),
                          yes = ifelse(test = m.sigs$log10OR_high > 0, yes = lim, no = -lim),
                          no = m.sigs$log10OR_high)

  m.sigs$log10OR_low = ifelse(test = is.infinite(m.sigs$log10OR_low),
                               yes = ifelse(test = m.sigs$log10OR_low > 0, yes = lim, no = -lim),
                               no = m.sigs$log10OR_low)

  if(!is.null(color)){
    color = color
    names(color) = c(m1Name, m2Name)
  }else{
    color = c("royalblue", "maroon")
    names(color) = c(m1Name, m2Name)
  }

  graphics::layout(mat = matrix(c(1, 2, 3, 4, 5, 5, 5, 5), byrow = TRUE, ncol = 4, nrow = 2), widths = c(4, 1, 1), heights = c(6, 1.2))
  par(mar = c(3, 1, 3, 5))
  plot(rep(0, nrow(m.sigs)), 1:nrow(m.sigs), xlim = c(-lim, lim), axes = FALSE, pch = NA, xlab = "", ylab = "", ylim = c(0.5, nrow(m.sigs)))
  segments(x0 = m.sigs$log10OR_low, y0 = 1:nrow(m.sigs), x1 = m.sigs$log10OR_high, y1 = 1:nrow(m.sigs), lwd = lineWidth, color[m.sigs$flow])
  points(m.sigs$log10OR, 1:nrow(m.sigs), pch = 16, cex = 0.5*(lineWidth))
  axis(side = 1, at = c(-lim, 0, lim), lwd = 2.2, font = 1, pos = 0.5, cex.axis = 1.3)
  segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = nrow(m.sigs)+0.2, col = "gray70", lwd = 2, lty = 2)
  mtext(text = m.sigs$Hugo_Symbol, side = 4, line = 0.2, at = 1:nrow(m.sigs),
        font = 3, las= 2, cex = geneFontSize, adj = 0)
  mtitle = paste(m2Name, ' (n = ', m2.sampleSize, ')', ' v/s ' , m1Name, ' (n = ' ,m1.sampleSize, ')', sep='')
  title(main = mtitle, font = 1, adj = 0, cex.main = titleSize)
  mtext(text = "Log odds ratio", side = 1, line = 2, font = 1, cex = 0.7*(titleSize))

  par(mar = c(3, 0, 3, 0))
  plot(rep(0, nrow(m.sigs)), 1:nrow(m.sigs), xlim = c(0, 1), axes = FALSE,
       pch = NA, xlab = "", ylab = "", ylim = c(0.5, nrow(m.sigs)))
  text(x = 0.5, y = 1:nrow(m.sigs), labels = as.numeric(unlist(m.sigs[,2])),
       adj = 0, font = 1, cex = 1.4*(geneFontSize))
  title(main = m1Name, cex.main = titleSize)

  par(mar = c(3, 0, 3, 0))
  plot(rep(0, nrow(m.sigs)), 1:nrow(m.sigs), xlim = c(0, 1), axes = FALSE,
       pch = NA, xlab = "", ylab = "", ylim = c(0.5, nrow(m.sigs)))
  text(x = 0.5, y = 1:nrow(m.sigs), labels = as.numeric(unlist(m.sigs[,3])),
       adj = 0, font = 1, cex = 1.4*(geneFontSize))
  title(main = m2Name, cex.main = titleSize)

  m.sigs$significance = ifelse(test =  as.numeric(m.sigs$pval) < 0.001, yes = "***", no =
                                 ifelse(test = as.numeric(m.sigs$pval) < 0.01, yes = "**", no =
                                          ifelse(test = as.numeric(m.sigs$pval) < 0.05, yes = "*", no = "NS")))

  par(mar = c(3, 0, 3, 0))
  plot(rep(0, nrow(m.sigs)), 1:nrow(m.sigs), xlim = c(0, 1), axes = FALSE,
       pch = NA, xlab = "", ylab = "", ylim = c(0.5, nrow(m.sigs)))
  text(x = 0.5, y = 1:nrow(m.sigs), labels = m.sigs$significance,
       adj = 0, font = 1, cex = 1.4*(geneFontSize))
  title(main = "p-value", cex.main = titleSize)


  plot.new()
  par(mar = c(0, 0, 0, 0))
  legend(x = "top", legend = names(color), lwd = lineWidth,
         col = color[c(m1Name, m2Name)], border = NA, bty = "n",
         cex = 1.1*(titleSize), horiz = TRUE, text.font = 1, xpd = TRUE, pch = 16)

}
