dashboard = function(maf, color = NULL, rmOutlier = TRUE, log_conv = FALSE, titv.color = NULL, sfs = statFontSize, fontSize = fs, n = 10, donut = pie, rawcount = TRUE, stat = NULL, titleSize = NULL, barcodes = NULL, barcodeSize = NULL){

  if(is.null(color)){
    #hard coded color scheme if user doesnt provide any
    col = get_vcColors()
  }else{
    col = color
  }

  vcs = getSampleSummary(maf)
  vcs = vcs[,colnames(vcs)[!colnames(x = vcs) %in% c('total', 'Amp', 'Del', 'CNV_total')], with = FALSE]
  vcs = vcs[,c(1,order(colSums(x = vcs[,2:(ncol(vcs)), with =FALSE]), decreasing = TRUE)+1), with =FALSE] #order based on most event
  vcs.m = data.table::melt(data = vcs, id = 'Tumor_Sample_Barcode')
  colnames(vcs.m) = c('Tumor_Sample_Barcode', 'Variant_Classification', 'N')

  data.table::setDF(vcs)
  rownames(x = vcs) = vcs$Tumor_Sample_Barcode
  vcs = vcs[,-1]
  vcs = t(vcs)

  lo = matrix(data = 1:6, nrow = 2, byrow = TRUE)
  layout(mat = lo, heights = c(3.5, 3), widths = c(3, 2, 2))
  par(cex.axis = fontSize, font = 3, cex.main = titleSize[1], lwd = 1.2)

  #--------------------------- variant classification plot -----------------
  vc.plot.dat = rev(rowSums(vcs))
  if(log_conv){
    vc.plot.dat = log10(vc.plot.dat)
  }

  xt = pretty(c(0, vc.plot.dat))

  par(mar = c(3, 9, 3, 1))
  b = barplot(vc.plot.dat, axes = FALSE, horiz = TRUE, col = col[names(vc.plot.dat)], border = NA,
              xlim = c(0, max(xt)), names.arg = rep("", length(vc.plot.dat)))
  abline(v = xt, h = 1:length(b)-0.25, lty = 2, lwd = 0.3, col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))
  axis(side = 2, at = b, labels = names(vc.plot.dat), lwd = 1.2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.9, font = 3, tick = FALSE)
  axis(side = 1, at = xt, lwd = 1.2, font = 3, las = 2, cex.axis = fontSize*0.9)
  title(main = "Variant Classification", adj = 0, cex.main = titleSize[1], font = 3)
  if(log_conv){
    axis(side = 2, at = 0, labels = "(log10)", lwd = 1.2, font = 3,
         las = 1, cex.axis = fontSize*0.9, hadj = 0.5, padj = 2, line = 0.75, tick = FALSE, outer = FALSE)
  }

  #--------------------------- variant type plot -----------------
  vt.plot.dat = maf@variant.type.summary
  vt.plot.dat = vt.plot.dat[,colnames(vt.plot.dat)[!colnames(x = vt.plot.dat) %in% c('total', 'CNV')], with = FALSE]
  vt.plot.dat = suppressWarnings(data.table::melt(vt.plot.dat[,c(2:(ncol(vt.plot.dat))), with = FALSE], id = NULL)[,sum(value), variable])
  colnames(vt.plot.dat)[2] = "sum"

  vt.cols = RColorBrewer::brewer.pal(n = 10, name = "Set3")
  if(log_conv){
    vt.plot.dat$sum = log10(vt.plot.dat$sum)
  }
  #xt = as.integer(seq(0, max(vt.plot.dat$sum), length.out = 4))
  xt = pretty(c(0, vt.plot.dat$sum))

  par(mar = c(3, 3, 3, 1))
  b = barplot(vt.plot.dat$sum, axes = FALSE, horiz = TRUE, col = vt.cols[1:length(vt.plot.dat$variable)],
              border = NA, xlim = c(0, max(xt)))
  abline(v = xt, h = 1:length(b)-0.25, lty = 2, lwd = 0.3, col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))
  axis(side = 2, at = b, labels = vt.plot.dat$variable, lwd = 1.2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.8, font = 3, tick = FALSE)
  axis(side = 1, at = xt, lwd = 1.2, font = 3, las = 2, cex.axis = fontSize*0.9)
  title(main = "Variant Type", adj = 0, cex.main = titleSize[1], font = 3)
  if(log_conv){
    axis(side = 2, at = 0, labels = "(log10)", lwd = 1.2, font = 3,
         las = 1, cex.axis = fontSize*0.9, hadj = 0.5, padj = 2, line = 0.75, tick = FALSE, outer = FALSE)
  }

  #--------------------------- titv summary plot -----------------
  titv = titv(maf = maf, useSyn = TRUE, plot = FALSE)
  titv.counts = titv$raw.counts
  titv.sums = data.table::melt(colSums(titv.counts[,2:7, with =FALSE]))
  titv.sums$class = rownames(titv.sums)
  if(!rawcount){
    titv.sums$raw_value = titv.sums$value
    titv.sums$value = titv.sums$value/sum(titv.sums$value)
    xt = seq(0, 1, 0.25)
  }else{
    xt = as.integer(seq(0, max(titv.sums$value, na.rm = TRUE), length.out = 4))
  }

  if(is.null(titv.color)){
    titv.color = get_titvCol()
  }

  par(mar = c(3, 3, 3, 1))
  b = barplot(titv.sums$value, axes = FALSE, horiz = TRUE, col = titv.color[rownames(titv.sums)],
              border = NA, xlim = c(0, xt[length(xt)]))
  axis(side = 2, at = b, labels = rownames(titv.sums), lwd = 1.2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.8, font = 3, tick = FALSE)
  axis(side = 1, at = xt, lwd = 1.2, font = 3, las = 2, cex.axis = fontSize*0.9)
  title(main = "SNV Class", adj = 0, cex.main = titleSize[1], font = 3)
  if(!rawcount){
    text(x = titv.sums$value+0.03, y = b, labels = titv.sums$raw_value,
         font = 4, col = "black", cex = fontSize, adj = 0)
  }
  abline(v = xt, h = 1:length(b)-0.25, lty = 2, lwd = 0.3, col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))

  #--------------------------- variant per sample plot -----------------

  if(barcodes){
    par(mar = c(6, 2, 3, 1))
  } else{
    par(mar = c(3, 2, 3, 1))
  }

  if(log_conv){
    vcs = apply(vcs, 2, function(x){
      x_fract = x / sum(x)
      x_log_total = log10(sum(x))
      x_fract * x_log_total
    })
    #Replace NaN's and remove those samples (resulting from samples with no variants but with CNVs)
    vcs[is.nan(x = vcs)] = 0
    vcs = vcs[,-which(colSums(x = vcs) == 0)]
  }

  b = barplot(vcs, col = col[rownames(vcs)], border = NA, axes = FALSE, names.arg =  rep("", ncol(vcs)))
  axis(side = 2, at = as.integer(seq(0, max(colSums(vcs)), length.out = 4)), lwd = 1.2, font = 3, las = 2,
       line = -0.3, hadj = 0.6, cex.axis = fontSize)
  title(main = "Variants per sample", adj = 0, cex.main = titleSize[1], font = 3, line = 2)

  if(barcodes){
    mtext(text = colnames(vcs), side = 1, line = 0.2, outer = FALSE, las = 2, at = b, cex = barcodeSize)
  }

  if(!is.null(stat)){
    if(stat == 'mean'){
      med.line = round(maf@summary[ID %in% "total", Mean], 2)
      if(log_conv){
        med.line = round(log10(med.line), 2)
      }
      df = data.frame(y = c(med.line), x = as.integer(0.8*nrow(getSampleSummary(maf))), label = c(paste('Mean: ', med.line, sep='')))
      title(main = paste0("Mean: ", med.line), adj = 0, cex.main = titleSize[1]*0.8, font = 3, line = 1)
      lines(x = c(1, b[length(b)]), y = c(med.line, med.line), col = "maroon", lwd = 1.2, lty = 2)
    }else if(stat == 'median'){
      med.line = round(maf@summary[ID %in% "total", Median], 2)
      if(log_conv){
        med.line = log10(med.line)
        title(main = paste0("Median: ", round(10^med.line)), adj = 0, cex.main = titleSize[1]*0.8, font = 3, line = 1)
      }else{
        title(main = paste0("Median: ", med.line), adj = 0, cex.main = titleSize[1]*0.8, font = 3, line = 1)
      }
      df = data.frame(y = c(med.line), x = as.integer(0.8*nrow(getSampleSummary(maf))), label = c(paste('Median: ', med.line, sep='')))
      lines(x = c(1, b[length(b)]), y = c(med.line, med.line), col = "maroon", lwd = 1.2, lty = 2)
    }
  }

  if(log_conv){
    mtext(text = "(log10)", side = 2, lwd = 1.2, font = 3, cex = fontSize*0.7, line = 1)
  }

  #--------------------------- vc summary plot -----------------
  par(mar = c(3, 2, 3, 1))
  boxH = vcs.m[,boxplot.stats(N)$stat[5], by = .(Variant_Classification)]
  colnames(boxH)[ncol(boxH)] = 'boxStat'
  bcol = col[levels(vcs.m$Variant_Classification)]
  #vcs.m$N = log10(vcs.m$N)
  b = boxplot(N ~ Variant_Classification, data = vcs.m, xaxt="n", outline=FALSE, lty=1, lwd = 1.4, outwex=0,
              staplewex=0, axes = FALSE, border = bcol)

  axis(side = 2, at = as.integer(seq(0, max(boxH[,boxStat], na.rm = TRUE), length.out = 4)),
       lwd = 1.2, font = 3, cex.axis = fontSize, las = 2)
  title(main = "Variant Classification \nsummary", adj = 0, cex.main = titleSize[1], font = 3, line = 1)
  abline(v = 1:length(bcol), h = as.integer(seq(0, max(boxH[,boxStat], na.rm = TRUE), length.out = 4)),
         lty = 2,
         lwd = 0.3, col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))

  #--------------------------- hugo-symbol plot -----------------
  gs = getGeneSummary(maf)
  nsamps = as.numeric(maf@summary[ID %in% "Samples", summary])
  gs.load = gs[,.(Hugo_Symbol, AlteredSamples)]
  gs.load[,AlteredSamples := round(AlteredSamples/nsamps, digits = 2) * 100]
  data.table::setDF(x = gs.load, rownames = gs.load$Hugo_Symbol)
  gs = gs[,colnames(gs)[!colnames(x = gs) %in% c('total', 'Amp', 'Del', 'CNV_total', 'MutatedSamples', 'AlteredSamples')], with = FALSE]

  if(nrow(gs) < n){
    gs.dat = gs
  }else{
    gs.dat = gs[1:n]
  }

  data.table::setDF(gs.dat)
  rownames(gs.dat) = gs.dat$Hugo_Symbol
  gs.dat = gs.dat[,-1]
  gs.dat = t(gs.dat)
  gs.dat = gs.dat[names(sort(rowSums(gs.dat), decreasing = TRUE)),, drop = FALSE]
  gs.dat = gs.dat[,names(sort(colSums(gs.dat))), drop = FALSE]

  xt = as.integer(seq(0, max(colSums(gs.dat))+2, length.out = 4))

  par(mar = c(3, 4, 3, 1))
  gs.load = gs.load[rev(colnames(gs.dat)),,]
  b = barplot(gs.dat, axes = FALSE, horiz = TRUE, col = col[rownames(gs.dat)], border = NA,
              xlim = c(0, max(xt)+(max(xt)*0.15)), names.arg = rep("", ncol(gs.dat)))
  axis(side = 2, at = b, labels = colnames(gs.dat), lwd = 1.2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.8, font = 3, tick = FALSE)
  axis(side = 1, at = xt, lwd = 1.2, font = 3, las = 2, cex.axis = fontSize*0.9)
  title(main = paste0('Top ',  n, '\nmutated genes'), adj = 0, cex.main = titleSize[1], font = 3)
  text(x = colSums(gs.dat)+1, y = b, labels = rev(paste0(gs.load$AlteredSamples, "%")),
       font = 4, col = "black", cex = fontSize*0.9, adj = 0)
  abline(h = b, v = xt,lty = 2, lwd = 0.3,
         col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))
}
