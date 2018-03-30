dashboard = function(maf, color = NULL, rmOutlier = TRUE, titv.color = NULL, sfs = statFontSize, fontSize = fs, n = 10, donut = pie, rawcount = TRUE, stat = NULL, titleSize = NULL){

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
  par(cex.axis = fontSize, font = 2, cex.main = titleSize[1], lwd = 2)

  #--------------------------- variant classification plot -----------------
  vc.plot.dat = rev(rowSums(vcs))
  xt = as.integer(seq(0, max(vc.plot.dat), length.out = 4))

  par(mar = c(3, 9, 3, 1))
  b = barplot(vc.plot.dat, axes = FALSE, horiz = TRUE, col = col[names(vc.plot.dat)], border = NA,
              xlim = c(0, max(xt)), names.arg = rep("", length(vc.plot.dat)))
  axis(side = 2, at = b, labels = names(vc.plot.dat), lwd = 2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.9, font = 2, tick = FALSE)
  axis(side = 1, at = xt, lwd = 2, font = 2, las = 2, cex.axis = fontSize*0.9)
  title(main = "Variant Classification", adj = 0, cex.main = titleSize[1], font = 2)

  #--------------------------- variant type plot -----------------
  vt.plot.dat = maf@variant.type.summary
  vt.plot.dat = vt.plot.dat[,colnames(vt.plot.dat)[!colnames(x = vt.plot.dat) %in% c('total', 'CNV')], with = FALSE]
  vt.plot.dat = suppressWarnings(data.table::melt(vt.plot.dat[,c(2:(ncol(vt.plot.dat))), with = FALSE], id = NULL)[,sum(value), variable])
  colnames(vt.plot.dat)[2] = "sum"

  vt.cols = RColorBrewer::brewer.pal(n = 10, name = "Set3")
  xt = as.integer(seq(0, max(vt.plot.dat$sum), length.out = 4))

  par(mar = c(3, 3, 3, 1))
  b = barplot(vt.plot.dat$sum, axes = FALSE, horiz = TRUE, col = vt.cols[1:length(vt.plot.dat$variable)],
              border = NA, xlim = c(0, max(xt)))
  axis(side = 2, at = b, labels = vt.plot.dat$variable, lwd = 2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.8, font = 2, tick = FALSE)
  axis(side = 1, at = xt, lwd = 2, font = 2, las = 2, cex.axis = fontSize*0.9)
  title(main = "Variant Type", adj = 0, cex.main = titleSize[1], font = 2)

  #--------------------------- titv summary plot -----------------
  titv = titv(maf = maf, useSyn = TRUE, plot = FALSE)
  titv.counts = titv$raw.counts
  titv.sums = data.table::melt(colSums(titv.counts[,2:7, with =FALSE]))
  titv.sums$class = rownames(titv.sums)
  if(!rawcount){
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
  axis(side = 2, at = b, labels = rownames(titv.sums), lwd = 2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.8, font = 2, tick = FALSE)
  axis(side = 1, at = xt, lwd = 2, font = 2, las = 2, cex.axis = fontSize*0.9)
  title(main = "SNV Class", adj = 0, cex.main = titleSize[1], font = 2)

  #--------------------------- variant per sample plot -----------------

  par(mar = c(3, 2, 3, 1))
  b = barplot(vcs, col = col[rownames(vcs)], border = NA, axes = FALSE, names.arg =  rep("", ncol(vcs)))
  axis(side = 2, at = as.integer(seq(0, max(colSums(vcs)), length.out = 4)), lwd = 2, font = 2, las = 2,
       line = -0.3, hadj = 0.6, cex.axis = fontSize)
  title(main = "Variants per sample", adj = 0, cex.main = titleSize[1], font = 2, line = 2)

  if(!is.null(stat)){
    if(stat == 'mean'){
      med.line = round(maf@summary[nrow(maf@summary),Mean], 2)
      df = data.frame(y = c(med.line), x = as.integer(0.8*nrow(getSampleSummary(maf))), label = c(paste('Mean: ', med.line, sep='')))
      title(main = paste0("Mean: ", med.line), adj = 0, cex.main = titleSize[1]*0.8, font = 2, line = 1)
    }else if(stat == 'median'){
      med.line = round(maf@summary[nrow(maf@summary),Median], 2)
      df = data.frame(y = c(med.line), x = as.integer(0.8*nrow(getSampleSummary(maf))), label = c(paste('Median: ', med.line, sep='')))
      title(main = paste0("Median: ", med.line), adj = 0, cex.main = titleSize[1]*0.8, font = 2, line = 1)
    }
  }else{
    med.line = round(max(maf@summary[,Median], na.rm = TRUE), 2)
    title(main = paste0("Median: ", med.line), adj = 0, cex.main = titleSize[1]*0.8, font = 2, line = 1)
  }


  lines(x = c(1, b[length(b)]), y = c(med.line, med.line), col = "maroon", lwd = 2, lty = 2)

  #--------------------------- vc summary plot -----------------
  par(mar = c(3, 2, 3, 1))
  boxH = vcs.m[,boxplot.stats(N)$stat[5], by = .(Variant_Classification)]
  colnames(boxH)[ncol(boxH)] = 'boxStat'
  bcol = col[levels(vcs.m$Variant_Classification)]
  b = boxplot(N ~ Variant_Classification, data = vcs.m, xaxt="n", outline=FALSE, lty=1, lwd = 1.4, outwex=0,
              staplewex=0, axes = FALSE, border = bcol)
#
#   boxplot(N ~ Variant_Classification, data = vcs.m,
#           col = col[levels(vcs.m$Variant_Classification)],
#           axes = FALSE, outline = FALSE, lwd = 1,
#           border = col[levels(vcs.m$Variant_Classification)],
#           staplewex=0, outwex=0)
  axis(side = 2, at = as.integer(seq(0, max(boxH[,boxStat], na.rm = TRUE), length.out = 4)),
       lwd = 2, font = 2, cex.axis = fontSize, las = 2)
  title(main = "Variant Classification \nsummary", adj = 0, cex.main = titleSize[1], font = 2, line = 1)

  #--------------------------- hugo-symbol plot -----------------
  gs = getGeneSummary(maf)
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

  xt = seq(0, max(colSums(gs.dat)), length.out = 4)

  par(mar = c(3, 4, 3, 1))
  b = barplot(gs.dat, axes = FALSE, horiz = TRUE, col = col[rownames(gs.dat)], border = NA, xlim = c(0, max(xt)), names.arg = rep("", ncol(gs.dat)))
  axis(side = 2, at = b, labels = colnames(gs.dat), lwd = 2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.8, font = 2, tick = FALSE)
  axis(side = 1, at = xt, lwd = 2, font = 2, las = 2, cex.axis = fontSize*0.9)
  title(main = paste0('Top ',  n, '\nmutated genes'), adj = 0, cex.main = titleSize[1], font = 2)

}

#   #--------------------------- oncoplot via ggplot -----------------
#
#   #Summary Table
#
#   hm.df = maf@oncoMatrix[1:10,]
#   hm.df = data.table::melt(hm.df)
#   #Reverse factors for decreasing order.
#   hm.df$Var1 = factor(hm.df$Var1, levels = rev(rownames(maf@oncoMatrix[1:20,])))
#
#   col = c(col, 'gray')
#   names(col)[length(col)] = ""
#
#   nsamples = as.numeric(maf.summary[ID %in% 'Samples',summary])
#
#   if(is.null(cohortName)){
#     cohortName = 'Cohort'
#   }
#
#   plot.title = paste(cohortName, " (n = ", nsamples,")", sep = '')
#
#   hm.gg = ggplot(data = hm.df, aes(x = Var2, y = Var1))+geom_tile(aes(fill = as.character(value)), color = 'white')+
#     scale_fill_manual(values = col)+
#     theme(legend.position = 'none', axis.line = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())+
#     xlab('')+ylab('')+ggtitle(label = plot.title)+theme(plot.background = element_rect(fill = '#F2F2F2'))
