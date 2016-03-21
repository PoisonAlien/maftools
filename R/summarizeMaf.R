summarizeMaf = function(maf){

  if('NCBI_Build' %in% colnames(maf)){
    NCBI_Build = unique(maf[,NCBI_Build])
    if(length(NCBI_Build) > 1){
      message('NOTE: Mutiple reference builds found!')
      NCBI_Build = do.call(paste, c(as.list(NCBI_Build), sep=";"))
      message(NCBI_Build)
    }
  }else{
    NCBI_Build = NA
  }

  if('Center' %in% colnames(maf)){
    Center = unique(maf[,Center])
    if(length(Center) > 1){
      message('Mutiple centers found.')
      Center = do.call(paste, c(as.list(Center), sep=";"))
      print(Center)
    }
  }else{
    Center = NA
  }

  #Variants per TSB
  tsb = maf[,.N, Tumor_Sample_Barcode]
  colnames(tsb)[2] = 'Variants'
  tsb = tsb[order(tsb$Variants, decreasing = T),]

  #summarise and casting by 'Variant_Classification'
  vc = maf[,.N, .(Tumor_Sample_Barcode, Variant_Classification )]
  vc.cast = dcast(data = vc, formula = Tumor_Sample_Barcode ~ Variant_Classification, fill = 0, value.var = 'N')
  vc.cast[,total:=rowSums(vc.cast[,2:ncol(vc.cast), with = F])]
  vc.cast = vc.cast[order(total, decreasing = T)]

  #summarise and casting by 'Variant_Type'
  vt = maf[,.N, .(Tumor_Sample_Barcode, Variant_Type )]
  vt.cast = dcast(data = vt, formula = Tumor_Sample_Barcode ~ Variant_Type, value.var = 'N', fill = 0)
  vt.cast[,total:=rowSums(vt.cast[,2:ncol(vt.cast), with = F])]
  vt.cast = vt.cast[order(total, decreasing = T)]

  #summarise and casting by 'Hugo_Symbol'
  hs = maf[,.N, .(Hugo_Symbol, Variant_Classification)]
  hs.cast = dcast(data = hs, formula = Hugo_Symbol ~Variant_Classification, fill = 0, value.var = 'N')
  hs.cast = hs.cast[order(rowSums(hs.cast[,2:ncol(hs.cast), with = F]), decreasing = T)] #ordering according to frequent mutated gene
  hs.cast[,total:=rowSums(hs.cast[,2:ncol(hs.cast), with = F])]
  #Get in how many samples a gene ismutated
  numMutatedSamples = maf[,.(MutatedSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
  #Merge and sort
  hs.cast = merge(hs.cast, numMutatedSamples, by = 'Hugo_Symbol')
  hs.cast = hs.cast[order(total, decreasing = T)]
  #Make a summarized table
  summary = data.table(ID = c('NCBI_Build', 'Center','Samples',colnames(vc.cast)[2:ncol(vc.cast)]), summary = c(NCBI_Build, Center,nrow(vc.cast), colSums(vc.cast[,2:ncol(vc.cast), with =F])))

  return(list(variants.per.sample = tsb, variant.type.summary = vt.cast, variant.classification.summary = vc.cast,
              gene.summary = hs.cast, summary = summary))
}
