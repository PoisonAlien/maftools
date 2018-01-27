dashboard = function(maf, color = NULL, rmOutlier = TRUE, titv.color = NULL, sfs = statFontSize, fontSize = fs, n = 10, donut = pie, rawcount = TRUE, stat = NULL){

  #--------------------------- Color code VC and gg theme -----------------

  db_theme = cowplot::theme_cowplot(font_size = fontSize, line_size = 0.75)+cowplot::background_grid(major = 'xy', minor = 'none')+
    theme(axis.text.x = element_text(face = 'bold', angle = 90, hjust = 1), axis.title.x = element_blank(), axis.text.y = element_text(face = 'bold'), axis.title.y = element_blank())+
    theme(legend.position = 'none')+
    theme(plot.title = element_text(color="black", face="bold", size = 9, hjust=0))+
    theme(plot.subtitle = element_text(color="#252525", face="bold", size = 8, hjust=0))



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


  if(!is.null(stat)){
    if(stat == 'mean'){
      med.line = round(maf@summary[nrow(maf@summary),Mean], 2)
      df = data.frame(y = c(med.line), x = as.integer(0.8*nrow(getSampleSummary(maf))), label = c(paste('Mean: ', med.line, sep='')))
    }else if(stat == 'median'){
      med.line = round(maf@summary[nrow(maf@summary),Median], 2)
      df = data.frame(y = c(med.line), x = as.integer(0.8*nrow(getSampleSummary(maf))), label = c(paste('Median: ', med.line, sep='')))
    }
  }else{
    med.line = round(max(maf@summary[,Median], na.rm = TRUE), 2)
  }

  vcs.gg = ggplot(data = vcs.m, aes(x = Tumor_Sample_Barcode, y = N, fill = Variant_Classification))+geom_bar(stat = 'identity')+scale_fill_manual(values = col)+
    ylab('# of Variants')+db_theme+theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())


  if(!is.null(stat)){
    vcs.gg = vcs.gg+geom_hline(yintercept = med.line, linetype = 3, size = 1, color = 'maroon')
    if(stat == 'mean'){
      vcs.gg = vcs.gg + ggtitle(label = 'Variants per sample', subtitle = paste0("Mean: ", med.line))+
        ggrepel::geom_label_repel(inherit.aes = FALSE, data = df, aes(x = x, y = y, label = label), fill = 'gray', fontface = 'bold', color = 'black', box.padding = unit(1, "lines"),
                                  point.padding = unit(1, "lines"), force = 10, size = sfs, nudge_y = 1)
    }else{
      vcs.gg = vcs.gg + ggtitle(label = 'Variants per sample', subtitle = paste0("Median: ", med.line))+
        ggrepel::geom_label_repel(inherit.aes = FALSE, data = df, aes(x = x, y = y, label = label), fill = 'gray', fontface = 'bold', color = 'black', box.padding = unit(1, "lines"),
                                  point.padding = unit(1, "lines"), force = 10, size = sfs, nudge_y = 1)
    }

  }else{
    vcs.gg = vcs.gg+ggtitle(label = 'Variants per sample', subtitle = paste0("Median: ", med.line))
  }

  #--------------------------- vc summary plot -----------------

  #bottom ggplot
  if(rmOutlier){
    #Get box heights from boxplot.srtas
    boxH = vcs.m[,boxplot.stats(N)$stat[5], by = .(Variant_Classification)]
    colnames(boxH)[ncol(boxH)] = 'boxStat'
    vcs.gg2 = ggplot(data = vcs.m, aes(x = Variant_Classification, y = N, fill = Variant_Classification))+
      geom_boxplot(outlier.shape = NA)+scale_fill_manual(values = col)+
      ylim(0, max(boxH$boxStat)+5)+db_theme

  } else{
    vcs.gg2 = ggplot(data = vcs.m, aes(x = Variant_Classification, y = N, fill = Variant_Classification))+
      geom_boxplot()+scale_fill_manual(values = col)+db_theme
  }

  vcs.gg2 = vcs.gg2+ggtitle(label = 'Variant Classification\nSummary')+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())


  #--------------------------- titv summary plot -----------------
  titv = titv(maf = maf, useSyn = TRUE, plot = FALSE)
  titv.counts = titv$raw.counts
  titv.sums = data.table::melt(colSums(titv.counts[,2:7, with =FALSE]))
  titv.sums$class = rownames(titv.sums)
  if(!rawcount){
    titv.sums$value = titv.sums$value/sum(titv.sums$value)
  }


  if(is.null(titv.color)){
    #titv.color = RColorBrewer::brewer.pal(n = 6, name = 'Set3')
    titv.color = c('coral4', 'lightcyan4', 'cornflowerblue', 'lightsalmon1', 'forestgreen', 'deeppink3')
    names(titv.color) = c('C>T', 'C>G', 'C>A', 'T>A', 'T>C', 'T>G')
  }


  titv.gg = ggplot(data = titv.sums, aes(x = class, y = value, fill = class))+
            geom_bar(stat = 'identity')+scale_fill_manual(values = titv.color)+
            coord_flip()+ggtitle(label = 'SNV Class')+db_theme

  if(!rawcount){
    titv.gg = titv.gg+scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), limits = c(0, 1))
  }



  #--------------------------- variant type plot -----------------
  vt.plot.dat = maf@variant.type.summary
  vt.plot.dat = vt.plot.dat[,colnames(vt.plot.dat)[!colnames(x = vt.plot.dat) %in% c('total', 'CNV')], with = FALSE]
  vt.plot.dat = suppressWarnings( data.table::melt(vt.plot.dat[,c(2:(ncol(vt.plot.dat))), with = FALSE], id = NULL))

  vt.gg = ggplot(data = vt.plot.dat, aes(x = variable, y = value, fill = variable))+
                geom_bar(stat = 'identity')+coord_flip()+
                ggtitle('Variant Type')+db_theme

  #--------------------------- variant classification plot -----------------
  vc.plot.dat = vcs[,2:ncol(vcs), with =FALSE]
  vc.lvl = sort(colSums(vc.plot.dat))
  vc.plot.dat = suppressWarnings( data.table::melt(vc.plot.dat))
  vc.plot.dat$variable = factor(x = vc.plot.dat$variable, levels = names(vc.lvl))

  vc.gg = ggplot(data = vc.plot.dat, aes(x = variable, y = value, fill = variable))+
      geom_bar(stat = 'identity')+scale_fill_manual(values = col)+coord_flip()+
      ggtitle('Variant Classification')+db_theme

#--------------------------- hugo-symbol plot -----------------
  gs = getGeneSummary(maf)
  gs = gs[,colnames(gs)[!colnames(x = gs) %in% c('total', 'Amp', 'Del', 'CNV_total', 'MutatedSamples', 'AlteredSamples')], with = FALSE]

  if(nrow(gs) < n){
    gs.dat = gs
  }else{
    gs.dat = gs[1:n]
  }

  gs.lvl = gs.dat[,Hugo_Symbol]
  gs.dat = suppressWarnings(data.table::melt(gs.dat))
  gs.dat$Hugo_Symbol = factor(x = gs.dat$Hugo_Symbol, levels = rev(gs.lvl))

  gs.gg = ggplot(data = gs.dat, aes(x = Hugo_Symbol, y = value, fill = variable))+
    geom_bar(stat = 'identity')+scale_fill_manual(values = col)+coord_flip()+
    ggtitle(label = paste0('Top ',  length(gs.lvl), '\nmutated genes'))+db_theme


  #--------------------------- Organize plots -----------------

  dash.gg = cowplot::plot_grid(vc.gg, vt.gg, titv.gg, vcs.gg, vcs.gg2, gs.gg,
                               nrow = 2, ncol = 3, rel_widths = c(1.2, 1, 1), rel_heights = c(1, 1.05)
                               )
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

