summarizeGistic = function(gistic){

  nGenes = length(unique(gistic[,Hugo_Symbol]))

  tsb = gistic[,.N, Tumor_Sample_Barcode]
  colnames(tsb)[2] = 'CNV'
  tsb = tsb[order(tsb$CNV, decreasing = TRUE),]

  vc = gistic[,.N, .(Tumor_Sample_Barcode, Variant_Classification )]
  vc.cast = data.table::dcast(data = vc, formula = Tumor_Sample_Barcode ~ Variant_Classification, fill = 0, value.var = 'N')

  vc.cast[,total:=rowSums(vc.cast[,2:ncol(vc.cast), with = FALSE])]
  vc.cast = vc.cast[order(total, decreasing = TRUE)]

  vc.mean = as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, mean))))
  vc.median = as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, median))))

  hs = gistic[,.N, .(Hugo_Symbol, Variant_Classification)]
  hs.cast = data.table::dcast(data = hs, formula = Hugo_Symbol ~Variant_Classification, fill = 0, value.var = 'N')
  hs.cast[,total:=rowSums(hs.cast[,2:ncol(hs.cast), with = FALSE])]
  hs.cast = hs.cast[order(total, decreasing = TRUE)]

  numAlteredSamples = gistic[,.(AlteredSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]

  hs.cast = merge(hs.cast, numAlteredSamples, by = 'Hugo_Symbol', all.x = TRUE)
  hs.cast = hs.cast[order(total, decreasing = TRUE)]

  summary = data.table::data.table(ID = c('Samples', 'nGenes', 'cytoBands', colnames(vc.cast)[2:ncol(vc.cast)]),
                                   summary = c(nrow(vc.cast), nGenes, length(gistic[,unique(as.character(Cytoband))]), colSums(vc.cast[,2:ncol(vc.cast), with =FALSE])))

  cytoBand = gistic[,.(nGenes = length(unique(Hugo_Symbol)), nSamples = length(unique(Tumor_Sample_Barcode)), Variant_Classification = unique(Variant_Classification)), by = Cytoband]

  return(list(summary = summary, cnv.summary = vc.cast,
              gene.summary = hs.cast, cytoband.summary = cytoBand))
}

gisticMap = function(gistic){
  dat = gistic

  # oncomat = data.table::dcast(data = dat, formula = Cytoband ~ Tumor_Sample_Barcode, value.var = 'Variant_Classification', fill = '' ,
  #       fun.aggregate = function(x) {ifelse(test = length(as.character(x))>1 ,no = as.character(x), yes = vcr(x))})

  oncomat = data.table::dcast(data = dat, formula = Cytoband ~ Tumor_Sample_Barcode, value.var = 'Variant_Classification', fill = '' ,
                              fun.aggregate = function(x) {
                                x = unique(as.character(x))
                                xad = x[x %in% c('Amp', 'Del')]
                                x = ifelse(test = length(xad) > 1, no = xad, yes = 'Complex')
                                return(x)
                              })

  #If maf contains only one sample converting to matrix is not trivial.
  if(ncol(oncomat) == 2){
    genes = oncomat[,Cytoband]
    sampleId = colnames(oncomat)[2]
    oncomat = as.matrix(data.frame(row.names = genes, sample = oncomat[,2, with =FALSE]))
  }else if(nrow(oncomat) == 1){
    #If MAF has only one gene
    gene = oncomat[,Cytoband]
    oncomat[,Cytoband:= NULL]
    oncomat = as.matrix(oncomat)
    rownames(oncomat) = gene
    sampleID = colnames(oncomat)
  }else{
    oncomat = as.matrix(oncomat)
    rownames(oncomat) = oncomat[,1]
    oncomat = oncomat[,-1]
  }

  variant.classes = as.character(unique(dat[,Variant_Classification]))
  variant.classes = c('',variant.classes)
  names(variant.classes) = 0:(length(variant.classes)-1)

  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes)){
    oncomat[oncomat == variant.classes[i]] = names(variant.classes)[i]
  }


  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = gene
    colnames(mdf) = sampleID
    return(list(oncomat = oncomat.copy, nummat = mdf, vc = variant.classes))
  }

  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)

  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId

    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId

    return(list(oncomat = oncomat.copy, nummat = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]

    mdf = mdf[, -ncol(mdf)]

    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix

    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort

    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy

    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]

    return(list(oncomat = oncomat.copy, nummat = mdf, vc = variant.classes))
  }
}
