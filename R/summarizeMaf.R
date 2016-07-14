
#Summarizing MAF
summarizeMaf = function(maf){

  maf.cnv = maf[Variant_Type %in% 'CNV']

  maf = maf[!Variant_Type %in% 'CNV']

  if('NCBI_Build' %in% colnames(maf)){
    NCBI_Build = unique(maf[,NCBI_Build])
    NCBI_Build = NCBI_Build[!is.na(NCBI_Build)]

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
    Center = Center[!is.na(Center)]
    if(length(Center) > 1){
      message('Mutiple centers found.')
      Center = do.call(paste, c(as.list(Center), sep=";"))
      print(Center)
    }
  }else{
    Center = NA
  }

  #nGenes
  nGenes = length(unique(maf[,Hugo_Symbol]))

  #Variants per TSB
  tsb = maf[,.N, Tumor_Sample_Barcode]
  colnames(tsb)[2] = 'Variants'
  tsb = tsb[order(tsb$Variants, decreasing = TRUE),]

  #summarise and casting by 'Variant_Classification'
  vc = maf[,.N, .(Tumor_Sample_Barcode, Variant_Classification )]
  vc.cast = data.table::dcast(data = vc, formula = Tumor_Sample_Barcode ~ Variant_Classification, fill = 0, value.var = 'N')
  vc.cast[,total:=rowSums(vc.cast[,2:ncol(vc.cast), with = FALSE])]
  vc.cast = vc.cast[order(total, decreasing = TRUE)]

  vc.mean = as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, mean))))
  vc.median = as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, median))))


  #summarise and casting by 'Variant_Type'
  vt = maf[,.N, .(Tumor_Sample_Barcode, Variant_Type )]
  vt.cast = data.table::dcast(data = vt, formula = Tumor_Sample_Barcode ~ Variant_Type, value.var = 'N', fill = 0)
  vt.cast[,total:=rowSums(vt.cast[,2:ncol(vt.cast), with = FALSE])]
  vt.cast = vt.cast[order(total, decreasing = TRUE)]

  #summarise and casting by 'Hugo_Symbol'
  hs = maf[,.N, .(Hugo_Symbol, Variant_Classification)]
  hs.cast = data.table::dcast(data = hs, formula = Hugo_Symbol ~Variant_Classification, fill = 0, value.var = 'N')
  hs.cast = hs.cast[order(rowSums(hs.cast[,2:ncol(hs.cast), with = FALSE]), decreasing = TRUE)] #ordering according to frequent mutated gene
  hs.cast[,total:=rowSums(hs.cast[,2:ncol(hs.cast), with = FALSE])]
  #Get in how many samples a gene ismutated
  numMutatedSamples = maf[,.(MutatedSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
  #Merge and sort
  hs.cast = merge(hs.cast, numMutatedSamples, by = 'Hugo_Symbol')
  hs.cast = hs.cast[order(MutatedSamples, decreasing = TRUE)]
  #Make a summarized table
  summary = data.table::data.table(ID = c('NCBI_Build', 'Center','Samples', 'nGenes',colnames(vc.cast)[2:ncol(vc.cast)]),
                       summary = c(NCBI_Build, Center, nrow(vc.cast), nGenes, colSums(vc.cast[,2:ncol(vc.cast), with =FALSE])))
  summary[,Mean := vc.mean]
  summary[,Median := vc.median]

  return(list(variants.per.sample = tsb, variant.type.summary = vt.cast, variant.classification.summary = vc.cast,
              gene.summary = hs.cast, summary = summary))
}

# This is using data.table. Very Fast :) :) :)
createOncoMatrix = function(maf){

    message('Creating oncomatrix (this might take a while)..')

     oncomat = data.table::dcast(data = maf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                                 fun.aggregate = function(x) {ifelse(test = length(as.character(x))>1 ,
                                no = as.character(x), yes = vcr(x, gis = FALSE))
                                 }, value.var = 'Variant_Classification', fill = '')

    #If maf contains only one sample converting to matrix is not trivial.
    if(ncol(oncomat) == 2){
      genes = oncomat[,Hugo_Symbol]
      sampleId = colnames(oncomat)[2]
      oncomat = as.matrix(data.frame(row.names = genes, sample = oncomat[,2, with =FALSE]))
    }else if(nrow(oncomat) == 1){
      #If MAF has only one gene
      gene = oncomat[,Hugo_Symbol]
      oncomat[,Hugo_Symbol:= NULL]
      oncomat = as.matrix(oncomat)
      rownames(oncomat) = gene
      sampleID = colnames(oncomat)
      }else{
      oncomat = as.matrix(oncomat)
      rownames(oncomat) = oncomat[,1]
      oncomat = oncomat[,-1]
      }

     variant.classes = as.character(unique(maf[,Variant_Classification]))
     variant.classes = c('',variant.classes, 'Multi_Hit')
     names(variant.classes) = 0:(length(variant.classes)-1)

     #Complex variant classes will be assigned a single integer.
     vc.onc = unique(unlist(apply(oncomat, 2, unique)))
     vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
     names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
     variant.classes2 = c(variant.classes, vc.onc)

     oncomat.copy <- oncomat
    #Make a numeric coded matrix
    for(i in 1:length(variant.classes2)){
      oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
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

    message('Sorting..')

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

#This is small function to sort genes according to total samples in which it is mutated.
sortByMutation = function(numMat, maf){

  geneOrder = getGeneSummary(x = maf)[,Hugo_Symbol]
  numMat = numMat[geneOrder[geneOrder %in% rownames(numMat)],]

  numMat[numMat != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
  tnumMat = t(numMat) #transposematrix
  numMat = t(tnumMat[do.call(order, c(as.list(as.data.frame(tnumMat)), decreasing = TRUE)), ]) #sort

  return(numMat)
}


# This is the original function for creating oncomatrix till 0.99.20. Verrrrrryyyyyy SLOWWWWWW !!!
#
# createOncoMatrix = function(maf){
#   message('Creating oncomatrix (this might take a while)..')
#
#   #tsbs
#   tsbs = levels(maf[,Tumor_Sample_Barcode])
#   #get all unique genes
#   genes = unique(maf[,Hugo_Symbol])
#
#   variant.classes = levels(maf[,Variant_Classification])
#   mdf = data.table(Hugo_Symbol = genes)
#
#   #prgress bar for for-loop
#   pb <- txtProgressBar(min = 0, max = length(tsbs), style = 3) #progress bar
#
#   for(i in 1:length(tsbs)) {
#     #message(paste('processing ',tsbs[i]))
#     barcode1 = maf[Tumor_Sample_Barcode == tsbs[i]]
#
#     if (nrow(barcode1) == 0)
#       break
#
#     x = barcode1[,.N, by = .(Hugo_Symbol, Variant_Classification)] #summerize by Hugo_Symbol and Variant_Classification
#
#     x.hits = x[,.N, Hugo_Symbol] #count number of mutation per gene
#     x.multihit.genes = x.hits[N > 1, Hugo_Symbol] #check for genes with >1 mutation
#     x.singlehits = x.hits[N == 1, Hugo_Symbol] #check for genes with single mutation
#
#     xs = x[Hugo_Symbol %in% x.singlehits, .(Hugo_Symbol, Variant_Classification)] #data table for single mutation
#
#     if(length(x.multihit.genes) > 0){
#       xm = data.table(Hugo_Symbol = x.multihit.genes, Variant_Classification = 'Multi_Hit') #data table for multi hit gene but Variant_Classification set to 'multi_hits'
#       x.hits.dt = rbind(xs, xm) #bind them togeather
#       colnames(x.hits.dt)[2] = tsbs[i] #change colnames to TSB
#       mdf = merge(mdf, x.hits.dt, by = 'Hugo_Symbol', all = TRUE) #merge
#     } else{
#       colnames(xs)[2] = tsbs[i]
#       mdf = merge(mdf, xs, by = 'Hugo_Symbol', all = TRUE) #merge
#     }
#     setTxtProgressBar(pb, i)
#   }
#
#   genes = as.character(mdf$Hugo_Symbol)
#   mdf = data.frame(mdf)
#   mdf = mdf[, -1]
#   mdf = data.frame(lapply(mdf, as.character), stringsAsFactors = FALSE)
#   mdf.copy = mdf #make a copy of this
#   mdf.copy[is.na(mdf.copy)] = ""
#   rownames(mdf.copy) = genes
#   mdf[is.na(mdf)] = 0
#
#   variant.classes = as.list(apply(mdf, 2, unique))
#   variant.classes = unique(as.vector(unlist(variant.classes)))
#   variant.classes = variant.classes[!variant.classes == "0"]
#
#   for (i in 1:length(variant.classes)) {
#     mdf[mdf == variant.classes[i]] = i
#   }
#
#   #numeric code for variant classes
#   names(variant.classes) = 1:length(variant.classes)
#
#   #convert dataframe to matrix
#   mdf = as.matrix(apply(mdf, 2, function(x) as.numeric(as.character(x))))
#   rownames(mdf) = genes #this is numeric coded matrix - maybe useful to used
#
#   message('Sorting..')
#   #Add total variants per gene
#   mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
#     length(x[x != "0"])
#   }))
#
#   #Sort the matrix
#   mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
#   colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
#   nMut = mdf[, ncol(mdf)]
#   mdf = mdf[, -ncol(mdf)]
#
#   mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix
#
#   mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
#   tmdf = t(mdf) #transposematrix
#   mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort
#
#   mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
#   mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
#   mdf = mdf.temp.copy
#
#   mat.numericCode = mdf #This is sorted, numeric coded matrix - maybe useful to user
#
#   mdf[is.na(mdf)] = 0
#   mdf.copy = mdf.copy[rownames(mdf), ]
#   #if names begin with numericasl they will be converted to ^X, sub it to get original id
#   colnames(mdf.copy) = gsub(pattern = '^X', replacement = '', x = colnames(mdf.copy))
#   mdf.copy = as.matrix(mdf.copy[, colnames(mdf)])
#
#   return(list(oncomat = mdf.copy, nummat = mat.numericCode, vc = variant.classes))
#
# }
