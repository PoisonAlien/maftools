mafCompare = function(m1, m2, m1Name = NULL, m2Name = NULL, minMut = 5, pVal = 0.05, showPval = FALSE){
 m1.gs = getGeneSummary(x = m1)
 m2.gs = getGeneSummary(x = m2)

 if(is.null(m1Name)){
   m1Name = 'M1'
 }

 if(is.null(m2Name)){
   m2Name = 'M2'
 }

 com.genes = intersect(m1.gs[,Hugo_Symbol], m2.gs[,Hugo_Symbol])

 m1.sampleSize = as.numeric(m1@summary[4,summary])
 m2.sampleSize = as.numeric(m2@summary[4,summary])

 m1.gs.comGenes = m1.gs[Hugo_Symbol %in% com.genes]
 m2.gs.comGenes = m2.gs[Hugo_Symbol %in% com.genes]

 m.gs.meged = merge(m1.gs.comGenes[,.(Hugo_Symbol, MutatedSamples)], m2.gs.comGenes[,.(Hugo_Symbol, MutatedSamples)],
                    by = 'Hugo_Symbol')
 m.gs.meged = as.data.frame(m.gs.meged)

 fisherTable = c()

 for(i in 1:nrow(m.gs.meged)){
   gene = m.gs.meged[i, 1]
   m1Mut = m.gs.meged[i,2]
   m2Mut = m.gs.meged[i,3]

   xf = fisher.test(matrix(c(m1Mut, m1.sampleSize, m2Mut, m2.sampleSize),
                           byrow = TRUE, nrow = 2), conf.int = TRUE, conf.level = 0.95)

   pval = xf$p.value
   or = xf$estimate
   ci.up = xf$conf.int[1]
   ci.low = xf$conf.int[2]
   tdat = data.table(Hugo_Symbol = gene, m1Mut , m2Mut, pval = pval, or = or, ci.up = ci.up, ci.low = ci.low)
   fisherTable = rbind(fisherTable, tdat)
 }

  fisherTable = fisherTable[order(pval)]

  colnames(fisherTable)[2:3] = c(m1Name, m2Name)

  m.sigs = fisherTable[pval < pVal]
  m.sigs$Hugo_Symbol = factor(x = m.sigs$Hugo_Symbol, levels = rev(m.sigs$Hugo_Symbol))
  m.sigs[,log10OR := log10(or)]
  m.sigs$label = paste('pval: ',round(m.sigs$pval, digits = 5), sep = '')
  m.sigs$flow = ifelse(test = m.sigs$log10OR < 0, yes = m2Name, no = m1Name)


  lim = max(abs(c(log10(m.sigs$ci.up), log10(m.sigs$ci.low))))

  gg.fp = ggplot(data = m.sigs, aes(x = Hugo_Symbol, y = log10OR, color = as.character(flow), label = label))+geom_point(size = 4)+
    geom_errorbar(aes(ymin = log10(ci.low), ymax = log10(ci.up)), size = 1, width = 0.20)+
    coord_flip()+cowplot::theme_cowplot(font_size = 9)+cowplot::background_grid(major = 'x')+
    geom_hline(yintercept = 0, linetype = 'dotted')+xlab('Gene')+ylab('log10 (Odds Ratio)')+
    theme(axis.line.y = element_blank(), legend.position = 'bottom', legend.title = element_blank())+ylim(-lim, lim)

  if(showPval){
    gg.fp+ggrepel::geom_label_repel(nudge_x = 0.2)
  }

  #gg.fp = gg.fp+annotate(geom = 'text', y = c(-2, 2), x = c(nrow(m.sigs)+2, nrow(m.sigs)+2), label = c(m1Name, m2Name))

 print(gg.fp)

  return(gg.fp)

}
