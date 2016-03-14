summarizeMaf = function(maf){
  #Variants per TSB
  tsb = maf[,.N, Tumor_Sample_Barcode]
  colnames(tsb)[2] = 'Variants'
  tsb = tsb[order(tsb$Variants, decreasing = T),]

  #summarise and casting by 'Variant_Classification'
  vc = maf[,.N, .(Tumor_Sample_Barcode, Variant_Classification )]
  vc.cast = dcast(data = vc, formula = Tumor_Sample_Barcode ~ Variant_Classification, fill = 0, value.var = 'N')
  vc.cast[,total:=rowSums(vc.cast[,2:ncol(vc.cast), with = F])]

  #summarise and casting by 'Variant_Type'
  vt = maf[,.N, .(Tumor_Sample_Barcode, Variant_Type )]
  vt.cast = dcast(data = vt, formula = Tumor_Sample_Barcode ~ Variant_Type, value.var = 'N', fill = 0)
  vt.cast[,total:=rowSums(vt.cast[,2:ncol(vt.cast), with = F])]

  #summarise and casting by 'Hugo_Symbol'
  hs = maf[,.N, .(Hugo_Symbol, Variant_Classification)]
  hs.cast = dcast(data = hs, formula = Hugo_Symbol ~Variant_Classification, fill = 0, value.var = 'N')
  hs.cast = hs.cast[order(rowSums(hs.cast[,2:ncol(hs.cast), with = F]), decreasing = T)] #ordering according to frequent mutated gene
  hs.cast[,total:=rowSums(hs.cast[,2:ncol(hs.cast), with = F])]

  #Make a summarized table
  summary = data.table(ID = c('Samples',colnames(vc.cast)[2:ncol(vc.cast)]), summary = c(nrow(vc.cast), colSums(vc.cast[,2:ncol(vc.cast), with =F])))

  return(list(variants.per.sample = tsb, variant.type.summary = vt.cast, variant.classification.summary = vc.cast,
              gene.summary = hs.cast, summary = summary))
}
