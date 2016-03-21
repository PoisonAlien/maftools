#' Read MAF files.
#' @description Takes tab delimited MAF (can be plain text or gz compressed) file as an input and summarizes it in various ways. Also creates oncomatrix - helpful for visualization.
#'
#' @param maf tab delimited MAF file. File can also be gz compressed.
#' @param removeSilent logical. Whether to discard silent (with no functional impact) mutations ("Silent","Intron","RNA","3'UTR"). Default is TRUE.
#' @param useAll logical. Whether to use all variants irrespective of values in Mutation_Status. Defaults to False. Only uses with values Somatic.
#' @return An object of class MAF.
#' @export


read.maf = function(maf, removeSilent = T, useAll = F){

  message('reading maf..')

  if(as.logical(length(grep(pattern = 'gz$', x = maf, fixed = F)))){
    #If system is Linux use fread, else use gz connection to read gz file.
    if(Sys.info()[['sysname']] == 'Linux'){
      maf = suppressWarnings(fread(input = paste('zcat <', maf), sep = '\t', stringsAsFactors = F, verbose = F, data.table = T, showProgress = T, header = T))
    } else{
      maf.gz = gzfile(description = maf, open = 'r')
      suppressWarnings(maf <- data.table(read.csv(file = maf.gz, header = T, sep = '\t', stringsAsFactors = F)))
      close(maf.gz)
    }
  } else{
    suppressWarnings(maf <- fread(input = maf, sep = "\t", stringsAsFactors = F, verbose = F, data.table = T, showProgress = T, header = T))
  }


  required.fields = c('Hugo_Symbol', 'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode') #necessary fields.
  missing.fileds = required.fields[!required.fields %in% colnames(maf)] #check if any of them is missing

  if(length(missing.fileds) > 0){
    stop(paste('missing', missing.fileds, 'from MAF\n')) #stop if any of required.fields are missing
  }

  #validation check for variants classified as Somatic in Mutation_Status field.
  if(length(colnames(maf)[colnames(x = maf) %in% 'Mutation_Status']) > 0){
    if(!useAll){
      message('Using only Somatic variants from Mutation_Status. Switch on useAll to include everything.')
      maf = subset(maf, Mutation_Status == 'Somatic')
    }else {
      message('Using all variants.')
    }
  }else{
    message('Mutation_Status not found. Assuming all variants are Somatic and validated.')
  }

  #convert "-" to "." in "Tumor_Sample_Barcode" to avoid complexity in naming
  maf$Tumor_Sample_Barcode = gsub(pattern = '-', replacement = '.', x = as.character(maf$Tumor_Sample_Barcode))

  #remove silent mutations with no/lesser functional impact
  silent = c("Silent", "Intron", "RNA", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR")
  maf.silent = maf[Variant_Classification %in% silent]

  if(removeSilent){

    if(nrow(maf.silent) > 0){
      maf.silent.vc = maf.silent[,.N, .(Tumor_Sample_Barcode, Variant_Classification)]
      maf.silent.vc.cast = dcast(data = maf.silent.vc, formula = Tumor_Sample_Barcode ~ Variant_Classification, fill = 0, value.var = 'N') #why dcast is not returning it as data.table ?
      summary.silent = data.table(ID = c('Samples',colnames(maf.silent.vc.cast)[2:ncol(maf.silent.vc.cast)]),
                                  N = c(nrow(maf.silent.vc.cast), colSums(maf.silent.vc.cast[,2:ncol(maf.silent.vc.cast), with = F])))

      maf = maf[!Variant_Classification %in% silent] #Remove silent variants from main table
      message(paste('removed',nrow(maf.silent), 'silent variants.'))
      print(summary.silent)
    } else{
      message(message(paste('removed',nrow(maf.silent), 'silent variants.')))
    }
  }else{
    message('Silent variants are being kept!')
  }

  #convert to factors
  maf$Variant_Type = as.factor(as.character(maf$Variant_Type))
  maf$Variant_Classification = as.factor(as.character(maf$Variant_Classification))
  maf$Tumor_Sample_Barcode = as.factor(as.character(maf$Tumor_Sample_Barcode))

  #tsbs
  tsbs = levels(maf[,Tumor_Sample_Barcode])
  #get all unique genes
  genes = unique(maf[,Hugo_Symbol])

  message('Summarizing..')
  mafSummary = summarizeMaf(maf = maf)
  print(mafSummary$summary)

  message("Frequently mutated genes..")
  print(mafSummary$gene.summary)

  message('Creating oncomatrix (this might take a while)..')
  variant.classes = levels(maf[,Variant_Classification])
  mdf = data.table(Hugo_Symbol = genes)

  #prgress bar for for-loop
  pb <- txtProgressBar(min = 0, max = length(tsbs), style = 3) #progress bar

  for(i in 1:length(tsbs)) {
    #message(paste('processing ',tsbs[i]))
    barcode1 = maf[Tumor_Sample_Barcode == tsbs[i]]

    if (nrow(barcode1) == 0)
      break

    x = barcode1[,.N, by = .(Hugo_Symbol, Variant_Classification)] #summerize by Hugo_Symbol and Variant_Classification

    x.hits = x[,.N, Hugo_Symbol] #count number of mutation per gene
    x.multihit.genes = x.hits[N > 1, Hugo_Symbol] #check for genes with >1 mutation
    x.singlehits = x.hits[N == 1, Hugo_Symbol] #check for genes with single mutation

    xs = x[Hugo_Symbol %in% x.singlehits, .(Hugo_Symbol, Variant_Classification)] #data table for single mutation

    if(length(x.multihit.genes) > 0){
      xm = data.table(Hugo_Symbol = x.multihit.genes, Variant_Classification = 'two_hit') #data table for multi hit gene but Variant_Classification set to 'two_hits'
      x.hits.dt = rbind(xs, xm) #bind them togeather
      colnames(x.hits.dt)[2] = tsbs[i] #change colnames to TSB
      mdf = merge(mdf, x.hits.dt, by = 'Hugo_Symbol', all = T) #merge
    } else{
      colnames(xs)[2] = tsbs[i]
      mdf = merge(mdf, xs, by = 'Hugo_Symbol', all = T) #merge
    }
    setTxtProgressBar(pb, i)
  }

  genes = as.character(mdf$Hugo_Symbol)
  mdf = data.frame(mdf)
  mdf = mdf[, -1]
  mdf = data.frame(lapply(mdf, as.character), stringsAsFactors = F)
  mdf.copy = mdf #make a copy of this
  mdf.copy[is.na(mdf.copy)] = ""
  rownames(mdf.copy) = genes
  mdf[is.na(mdf)] = 0

  variant.classes = as.list(apply(mdf, 2, unique))
  variant.classes = unique(as.vector(unlist(variant.classes)))
  variant.classes = variant.classes[!variant.classes == "0"]

  for (i in 1:length(variant.classes)) {
    mdf[mdf == variant.classes[i]] = i
  }

  #numeric code for variant classes
  names(variant.classes) = 1:length(variant.classes)

  #convert dataframe to matrix
  mdf = as.matrix(apply(mdf, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = genes #this is numeric coded matrix - maybe useful to used

  message('Sorting..')
  #Add total variants per gene
  mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
    length(x[x != "0"])
  }))

  #Sort the matrix
  mdf = mdf[order(mdf[, ncol(mdf)], decreasing = T), ]
  colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
  nMut = mdf[, ncol(mdf)]
  mdf = mdf[, -ncol(mdf)]

  mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix

  mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
  tmdf = t(mdf) #transposematrix
  mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = T)), ]) #sort

  mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
  mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
  mdf = mdf.temp.copy

  mat.numericCode = mdf #This is sorted, numeric coded matrix - maybe useful to user

  mdf[is.na(mdf)] = 0
  mdf.copy = mdf.copy[rownames(mdf), ]
  #if names begin with numericasl they will be converted to ^X, sub it to get original id
  colnames(mdf.copy) = gsub(pattern = '^X', replacement = '', x = colnames(mdf.copy))
  mdf.copy = as.matrix(mdf.copy[, colnames(mdf)])
  message('Done !')


  m = MAF(data = maf, variants.per.sample = mafSummary$variants.per.sample, variant.type.summary = mafSummary$variant.type.summary,
              variant.classification.summary = mafSummary$variant.classification.summary,gene.summary = mafSummary$gene.summary,
              oncoMatrix = mdf.copy, numericMatrix = mat.numericCode,summary = mafSummary$summary,
              classCode = variant.classes, maf.silent = maf.silent)

  return(m)

}
