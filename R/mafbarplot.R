#' Creates a bar plot
#' @description Takes an MAF object and generates a barplot of mutated genes color coded for variant classification
#' @param maf an \code{\link{MAF}} object
#' @param n Number of genes to include. Default 20.
#' @param genes Manually provide names of genes. Default NULL.
#' @param color named vector of colors for each Variant_Classification. Default NULL.
#' @param fontSize Default 0.7
#' @param includeCN Include copy number events if available? Default FALSE
#' @param legendfontSize Default 0.7
#' @param borderCol Default "#34495e". Set to `NA` for no border color.
#' @param showPct Default TRUE. Show percent altered samples.
#' @export
#' @examples
#' laml.maf = system.file("extdata", "tcga_laml.maf.gz", package = "maftools") #MAF file
#' laml = read.maf(maf = laml.maf)
#' mafbarplot(maf = laml)
#'
mafbarplot = function(maf, n = 20, genes = NULL, color = NULL, fontSize = 0.7, includeCN = FALSE, legendfontSize = 0.7, borderCol = "#34495e", showPct = TRUE){

  if(is.null(color)){
    #hard coded color scheme if user doesnt provide any
    col = get_vcColors()
  }else{
    col = color
  }

  gs = getGeneSummary(maf)
  nsamps = as.numeric(maf@summary[ID %in% "Samples", summary])
  gs.load = gs[,.(Hugo_Symbol, AlteredSamples)]
  gs.load[,AlteredSamples := round(AlteredSamples/nsamps, digits = 2) * 100]
  data.table::setDF(x = gs.load, rownames = gs.load$Hugo_Symbol)
  if(includeCN){
    gs = gs[,colnames(gs)[!colnames(x = gs) %in% c('total', 'CNV_total', 'MutatedSamples', 'AlteredSamples')], with = FALSE]
  }else{
    gs = gs[,colnames(gs)[!colnames(x = gs) %in% c('total', 'Amp', 'Del', 'CNV_total', 'MutatedSamples', 'AlteredSamples')], with = FALSE]
  }

  if(is.null(genes)){
    if(nrow(gs) < n){
      gs.dat = gs
    }else{
      gs.dat = gs[1:n]
    }
  }else{
    gs.dat = gs[Hugo_Symbol %in% genes]
    if(nrow(gs.dat) == 0){
      stop("Selected genes were not found!")
    }
  }

  data.table::setDF(gs.dat)
  rownames(gs.dat) = gs.dat$Hugo_Symbol
  gs.dat = gs.dat[,-1, drop = FALSE]

  gs.dat = t(gs.dat)
  gs.dat = gs.dat[names(sort(rowSums(gs.dat), decreasing = TRUE)),, drop = FALSE]
  gs.dat = gs.dat[,names(sort(colSums(gs.dat))), drop = FALSE]

  xt = as.integer(seq(0, max(colSums(gs.dat))+2, length.out = 4))

  par(mar = c(3, 5, 0.5, 1))
  gs.load = gs.load[rev(colnames(gs.dat)),,]
  b = barplot(gs.dat, axes = FALSE, horiz = TRUE, col = col[rownames(gs.dat)], border = borderCol,
              xlim = c(0, max(xt)+(max(xt)*0.15)), names.arg = rep("", ncol(gs.dat)))
  axis(side = 2, at = b, labels = colnames(gs.dat), lwd = 1.2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.8, font = 3, tick = FALSE)
  axis(side = 1, at = xt, lwd = 1.2, font = 1, las = 1, cex.axis = 1)
  #title(main = paste0('Top ',  n, '\nmutated genes'), adj = 0, cex.main = titleSize[1], font = 3)
  if(showPct){
    text(x = colSums(gs.dat), y = b, labels = rev(paste0(gs.load$AlteredSamples, "%")),
         font = 3, col = "black", cex = fontSize*0.9, adj = 0, xpd = TRUE, pos = 4)
  }

  abline(h = b, v = xt,lty = 2, lwd = 0.3,
         col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))

  leg_cols = col[rownames(gs.dat)]
  mtext(text = "No. of alterations", side = 1, line = 2, adj = 0)
  legend(x = "bottomright", legend = names(leg_cols), col = leg_cols, pch = 19, ncol = 1, cex = legendfontSize)
}
