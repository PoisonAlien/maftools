dashboard = function(maf, color = NULL, rmOutlier = TRUE, titv.color = NULL, sfs = statFontSize, fontSize = fs, n = 10){

  #--------------------------- Color code for VC -----------------


  if(is.null(color)){
    #hard coded color scheme if user doesnt provide any
    col = c(RColorBrewer::brewer.pal(12,name = "Paired"), RColorBrewer::brewer.pal(11,name = "Spectral")[1:3],'black')
    names(col) = names = c('Nonstop_Mutation','Frame_Shift_Del','IGR','Missense_Mutation','Silent','Nonsense_Mutation',
                           'RNA','Splice_Site','Intron','Frame_Shift_Ins','Nonstop_Mutation','In_Frame_Del','ITD','In_Frame_Ins','Translation_Start_Site',"Multi_Hit")
  }else{
    col = color
  }

  #--------------------------- variant per sample plot -----------------

  vcs = getSampleSummary(maf)
  vcs = vcs[,colnames(vcs)[!colnames(x = vcs) %in% c('total', 'Amp', 'Del', 'CNV_total')], with = FALSE]

  #vcs = vcs[,2:(ncol(vcs)-1), with = F]
  #vcs = vcs[order(rowSums(vcs[,2:ncol(vcs), with = F]), decreasing = T)] #order tsbs based on number of mutations
  vcs = vcs[,c(1,order(colSums(x = vcs[,2:(ncol(vcs)), with =FALSE]), decreasing = TRUE)+1), with =FALSE] #order based on most event

  #melt data frame
  vcs.m = data.table::melt(data = vcs, id = 'Tumor_Sample_Barcode')
  colnames(vcs.m) = c('Tumor_Sample_Barcode', 'Variant_Classification', 'N')

  vcs.m$Tumor_Sample_Barcode = factor(vcs.m$Tumor_Sample_Barcode,levels = vcs$Tumor_Sample_Barcode) #reorder tsb levels
  vcs.m$Variant_Classification = factor(x = vcs.m$Variant_Classification, levels = colnames(vcs)) #reorder Variant classification

  #For now keep it at 90
  textAngle = 90

  med.line = round(max(maf@summary[,Median], na.rm = TRUE), 2)
  df = data.frame(y = c(med.line), x = as.integer(0.8*nrow(getSampleSummary(maf))), label = c(paste('Median: ', med.line, sep='')))

  vcs.gg = ggplot(data = vcs.m, aes(x = Tumor_Sample_Barcode, y = N, fill = Variant_Classification))+geom_bar(stat = 'identity')+scale_fill_manual(values = col)+cowplot::theme_cowplot(font_size = fontSize)+
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), legend.position = 'none')+
    ylab('# of Variants')+xlab('Samples')+
    cowplot::background_grid(major = 'xy', minor = 'none')+geom_hline(yintercept = med.line, linetype = 3)+
    ggrepel::geom_label_repel(inherit.aes = FALSE, data = df, aes(x = x, y = y, label = label), fill = 'gray', fontface = 'bold', color = 'black', box.padding = unit(1, "lines"),
                              point.padding = unit(1, "lines"), force = 20, size = sfs, nudge_y = 2)

  vcs.gg = vcs.gg+ggtitle(label = 'Variants per sample')

  #--------------------------- vc summary plot -----------------

  #bottom ggplot
  if(rmOutlier){
    #Get box heights from boxplot.srtas
    boxH = vcs.m[,boxplot.stats(N)$stat[5], by = .(Variant_Classification)]
    colnames(boxH)[ncol(boxH)] = 'boxStat'
    vcs.gg2 = ggplot(data = vcs.m, aes(x = Variant_Classification, y = N, fill = Variant_Classification))+geom_boxplot(outlier.shape = NA)+scale_fill_manual(values = col)+
      cowplot::theme_cowplot(font_size = fontSize)+theme(axis.ticks = element_blank(), axis.text.x = element_blank(), legend.position = 'none')+
      theme(legend.key.size = unit(0.5, "cm"), legend.title = element_blank(), legend.text = element_text(size = 8))+ylim(0, max(boxH$boxStat)+5)+cowplot::background_grid(major = 'xy', minor = 'none')
  } else{
    vcs.gg2 = ggplot(data = vcs.m, aes(x = Variant_Classification, y = N, fill = Variant_Classification))+geom_boxplot()+scale_fill_manual(values = col)+cowplot::theme_cowplot(font_size = fontSize)+
      theme(axis.ticks = element_blank(), axis.text.x = element_blank(), legend.position = 'none')+
      theme(legend.key.size = unit(0.5, "cm"), legend.title = element_blank(), legend.text = element_text(size = 8))+cowplot::background_grid(major = 'xy', minor = 'none')
  }

  vcs.gg2 = vcs.gg2+ggtitle('Variant Classification Summary')


  #--------------------------- titv summary plot -----------------
  titv = titv(maf = maf, useSyn = TRUE, plot = FALSE)
  titv.counts = titv$raw.counts
  titv.sums = data.table::melt(colSums(titv.counts[,2:7, with =FALSE]))
  titv.sums$class = rownames(titv.sums)
  titv.sums$class = gsub(pattern = '-', replacement = '>', x = titv.sums$class)


  if(is.null(titv.color)){
    titv.color = RColorBrewer::brewer.pal(n = 6, name = 'Set3')
    names(titv.color) = c('C>T', 'C>G', 'C>A', 'T>A', 'T>C', 'T>G')
  }

  titv.gg = ggplot(data = titv.sums, aes(x = class, y = value, fill = class))+geom_bar(stat = 'identity')+scale_fill_manual(values = titv.color)+
    cowplot::theme_cowplot(font_size = fontSize)+cowplot::background_grid(major = 'xy')+coord_flip()+
    theme(legend.position = 'none')+xlab('')+ylab('N')+ggtitle('SNV Class')


  #--------------------------- variant type plot -----------------
  vt.plot.dat = maf@variant.type.summary
  vt.plot.dat = vt.plot.dat[,colnames(vt.plot.dat)[!colnames(x = vt.plot.dat) %in% c('total', 'CNV')], with = FALSE]
  vt.plot.dat = suppressWarnings( data.table::melt(vt.plot.dat[,c(2:(ncol(vt.plot.dat))), with = FALSE], id = NULL))

  vt.gg = ggplot(data = vt.plot.dat, aes(x = variable, y = value, fill = variable))+
    geom_bar(stat = 'identity')+coord_flip()+cowplot::theme_cowplot(font_size = fontSize)+
    cowplot::background_grid(major = 'xy')+
    theme(legend.position = 'none')+xlab('')+ylab('N')+ggtitle('Variant Type')

  #--------------------------- variant classification plot -----------------
#   vc.plot.dat = vcs
#   vc.plot.dat = vc.plot.dat[,c(2:(ncol(vc.plot.dat))), with =FALSE]
  vc.plot.dat = vcs[,2:ncol(vcs), with =FALSE]
  vc.lvl = sort(colSums(vc.plot.dat))
  vc.plot.dat = suppressWarnings( data.table::melt(vc.plot.dat))
  vc.plot.dat$variable = factor(x = vc.plot.dat$variable, levels = names(vc.lvl))

  vc.gg = ggplot(data = vc.plot.dat, aes(x = variable, y = value, fill = variable))+
    geom_bar(stat = 'identity')+scale_fill_manual(values = col)+coord_flip()+
    cowplot::theme_cowplot(font_size = fontSize)+cowplot::background_grid(major = 'xy')+
    theme(legend.position = 'none')+xlab('')+ylab('N')

  vc.gg = vc.gg+ggtitle('Variant Classification')


#--------------------------- hugo-symbol plot -----------------
  gs = getGeneSummary(maf)
  gs = gs[,colnames(gs)[!colnames(x = gs) %in% c('total', 'Amp', 'Del', 'CNV_total', 'MutatedSamples')], with = FALSE]
  gs.dat = gs[1:n]
  gs.lvl = gs.dat[,Hugo_Symbol]
  gs.dat = suppressWarnings(data.table::melt(gs.dat))
  gs.dat$Hugo_Symbol = factor(x = gs.dat$Hugo_Symbol, levels = rev(gs.lvl))

  gs.gg = ggplot(data = gs.dat, aes(x = Hugo_Symbol, y = value, fill = variable))+geom_bar(stat = 'identity')+scale_fill_manual(values = col)+coord_flip()+
    cowplot::theme_cowplot(font_size = fontSize)+cowplot::background_grid(major = 'xy')+
    theme(legend.position = 'none')+xlab('')+ylab('N')
  gs.gg = gs.gg+ggtitle('Frequently Mutated Genes')


  #--------------------------- Organize plots -----------------

  dash.gg = cowplot::plot_grid(vt.gg, vc.gg, titv.gg, vcs.gg, vcs.gg2, gs.gg)
  return(dash.gg)

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

