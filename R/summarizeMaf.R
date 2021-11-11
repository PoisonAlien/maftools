summarizeMaf = function(maf, anno = NULL, chatty = TRUE){

  if(!is.data.frame(maf)){
    #try to coerce into data.frame
    maf = data.table::as.data.table(maf)
  }

  if('NCBI_Build' %in% colnames(maf)){
    NCBI_Build = unique(maf[!Variant_Type %in% 'CNV', NCBI_Build])
    NCBI_Build = NCBI_Build[!is.na(NCBI_Build)]

    if(chatty){
      if(length(NCBI_Build) > 1){
        cat('--Mutiple reference builds found\n')
        NCBI_Build = do.call(paste, c(as.list(NCBI_Build), sep=";"))
        cat(NCBI_Build)
      }
    }
  }else{
    NCBI_Build = NA
  }

  if('Center' %in% colnames(maf)){
    Center = unique(maf[!Variant_Type %in% 'CNV', Center])
    #Center = Center[is.na(Center)]
    if(length(Center) > 1){
      Center = do.call(paste, c(as.list(Center), sep=";"))
      if(chatty){
        cat('--Mutiple centers found\n')
        cat(Center)
      }
    }
  }else{
    Center = NA
  }


  #nGenes
  nGenes = length(unique(maf[,Hugo_Symbol]))
  maf.tsbs = levels(maf[,Tumor_Sample_Barcode])
  nSamples = length(levels(maf$Tumor_Sample_Barcode))

  #Top 20 FLAGS - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267152/
  flags = flags(top = 20)

  #Variants per TSB
  tsb = maf[,.N, Tumor_Sample_Barcode]
  colnames(tsb)[2] = 'Variants'
  tsb = tsb[order(tsb$Variants, decreasing = TRUE),]

  #summarise and casting by 'Variant_Classification'
  vc = maf[,.N, .(Tumor_Sample_Barcode, Variant_Classification )]
  vc.cast = data.table::dcast(data = vc, formula = Tumor_Sample_Barcode ~ Variant_Classification, fill = 0, value.var = 'N', drop = FALSE)

  if(any(colnames(vc.cast) %in% c('Amp', 'Del'))){
    vc.cast.cnv = vc.cast[,c('Tumor_Sample_Barcode', colnames(vc.cast)[colnames(vc.cast) %in% c('Amp', 'Del')]), with =FALSE]
    vc.cast.cnv$CNV_total = rowSums(vc.cast.cnv[,2:ncol(vc.cast.cnv)], na.rm = TRUE)

    vc.cast = vc.cast[,!colnames(vc.cast)[colnames(vc.cast) %in% c('Amp', 'Del')], with =FALSE]
    vc.cast[,total:=rowSums(vc.cast[,2:ncol(vc.cast), with = FALSE])]

    vc.cast = merge(vc.cast, vc.cast.cnv, by = 'Tumor_Sample_Barcode', all = TRUE)[order(total, CNV_total, decreasing = TRUE)]

    vc.mean = as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, mean))))
    vc.median = as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, median))))

  }else{
    vc.cast = vc.cast[,total:=rowSums(vc.cast[,2:ncol(vc.cast), with = FALSE])][order(total, decreasing = TRUE)]

    vc.mean = round(as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, mean)))), 3)
    vc.median = round(as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, median)))), 3)
  }

  #summarise and casting by 'Variant_Type'
  vt = maf[,.N, .(Tumor_Sample_Barcode, Variant_Type )]
  vt.cast = data.table::dcast(data = vt, formula = Tumor_Sample_Barcode ~ Variant_Type, value.var = 'N', fill = 0, drop = FALSE)

  if(any(colnames(vt.cast) %in% c('CNV'))){
    vt.cast.cnv = vt.cast[,c('Tumor_Sample_Barcode', colnames(vt.cast)[colnames(vt.cast) %in% c('CNV')]), with =FALSE]

    vt.cast = vt.cast[,!colnames(vt.cast)[colnames(vt.cast) %in% c('CNV')], with =FALSE]
    vt.cast = vt.cast[,total:=rowSums(vt.cast[,2:ncol(vt.cast), with = FALSE])]

    vt.cast = merge(vt.cast, vt.cast.cnv, by = 'Tumor_Sample_Barcode', all = TRUE)[order(total, CNV, decreasing = TRUE)]
  }else{
    vt.cast = vt.cast[,total:=rowSums(vt.cast[,2:ncol(vt.cast), with = FALSE])][order(total, decreasing = TRUE)]
  }

  #summarise and casting by 'Hugo_Symbol'
  hs = maf[,.N, .(Hugo_Symbol, Variant_Classification)]
  hs.cast = data.table::dcast(data = hs, formula = Hugo_Symbol ~Variant_Classification, fill = 0, value.var = 'N')
  #----
  if(any(colnames(hs.cast) %in% c('Amp', 'Del'))){
    hs.cast.cnv = hs.cast[,c('Hugo_Symbol', colnames(hs.cast)[colnames(hs.cast) %in% c('Amp', 'Del')]), with = FALSE]
    hs.cast.cnv$CNV_total = rowSums(x = hs.cast.cnv[,2:ncol(hs.cast.cnv), with = FALSE], na.rm = TRUE)

    hs.cast = hs.cast[,!colnames(hs.cast)[colnames(hs.cast) %in% c('Amp', 'Del')], with = FALSE]
    hs.cast[,total:=rowSums(hs.cast[,2:ncol(hs.cast), with = FALSE], na.rm = TRUE)]

    hs.cast = merge(hs.cast, hs.cast.cnv, by = 'Hugo_Symbol', all = TRUE)[order(total, CNV_total, decreasing = TRUE)]
  }else{
    hs.cast[,total:=rowSums(hs.cast[,2:ncol(hs.cast), with = FALSE])]
    hs.cast = hs.cast[order(total, decreasing = TRUE)]
  }
  #----

  #Get in how many samples a gene ismutated
  numMutatedSamples = maf[!Variant_Type %in% 'CNV', .(MutatedSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
  numAlteredSamples = maf[, .(AlteredSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
  numAlteredSamples = merge(numMutatedSamples, numAlteredSamples, by = 'Hugo_Symbol', all = TRUE)
  #Merge and sort
  hs.cast = merge(hs.cast, numAlteredSamples, by = 'Hugo_Symbol', all = TRUE)[order(MutatedSamples, total, decreasing = TRUE)]
  #Replace NAs with 0
  hs.cast$AlteredSamples = ifelse(test = is.na(x = hs.cast$AlteredSamples), yes = 0, no = hs.cast$AlteredSamples)
  hs.cast$MutatedSamples = ifelse(test = is.na(x = hs.cast$MutatedSamples), yes = 0, no = hs.cast$MutatedSamples)
  #Make a summarized table
  summary = data.table::data.table(ID = c('NCBI_Build', 'Center','Samples', 'nGenes',colnames(vc.cast)[2:ncol(vc.cast)]),
                       summary = c(NCBI_Build, Center, nSamples, nGenes, colSums(vc.cast[,2:ncol(vc.cast), with =FALSE])))
  summary[,Mean := vc.mean]
  summary[,Median := vc.median]

  # if(chatty){
  #   print(summary)
  #
  #   cat("Gene Summary:\n")
  #   print(hs.cast)
  # }

  #Check for flags.
  if(nrow(hs.cast) > 10){
    topten = as.character(hs.cast[1:10, Hugo_Symbol])
    topten = topten[topten %in% flags]
      if(chatty){
        if(length(topten) > 0){
          cat('--Possible FLAGS among top ten genes:\n')
          for(temp in topten){
            cat(paste0("  ", temp, "\n"))
          }
      }
    }
  }

  if(chatty){
    cat("-Processing clinical data\n")
  }

  if(is.null(anno)){
    if(chatty){
      cat("--Missing clinical data\n")
    }
    sample.anno = data.table::data.table(Tumor_Sample_Barcode = maf.tsbs)
  }else if(is.data.frame(x = anno)){
    sample.anno = data.table::copy(x = anno)
    data.table::setDT(sample.anno)
      if(!'Tumor_Sample_Barcode' %in% colnames(sample.anno)){
        message(paste0('Available fields in provided annotations..'))
        print(colnames(sample.anno))
        stop(paste0('Tumor_Sample_Barcode column not found in provided clinical data. Rename column containing sample names to Tumor_Sample_Barcode if necessary.'))
      }
    }else{
      if(file.exists(anno)){
        sample.anno = data.table::fread(anno, stringsAsFactors = FALSE)
        if(!'Tumor_Sample_Barcode' %in% colnames(sample.anno)){
          message(paste0('Available fields in ', basename(anno), '..'))
          print(colnames(sample.anno))
          stop(paste0('Tumor_Sample_Barcode column not found in provided clinical data. Rename column name containing sample names to Tumor_Sample_Barcode if necessary.'))
        }
      }
    }

  #clean up annotation data
  colnames(sample.anno) = gsub(pattern = ' ', replacement = '_', x = colnames(sample.anno), fixed = TRUE) #replace spaces in column names for annotation data
  if(nrow(sample.anno) == 1){
    temp_colnames = colnames(sample.anno)
    sample.anno = as.data.frame(apply(sample.anno, 2, function(y) trimws(y))) #remove trailing whitespaces
    sample.anno = data.frame(t(unlist(sample.anno, use.names = FALSE)))
    colnames(sample.anno) = temp_colnames
  }else{
    sample.anno = as.data.frame(apply(sample.anno, 2, function(y) trimws(y))) #remove trailing whitespaces
  }

  sample.anno[sample.anno == ""] = NA #Replace blanks with NA
  #sample.anno = as.data.frame(apply(sample.anno, 2, function(y) gsub(pattern = " ", replacement = "_", x = y))) #replace spaces with _
  data.table::setDT(x = sample.anno)
  if(ncol(sample.anno) == 1){
    colnames(sample.anno)[1] = c("Tumor_Sample_Barcode")
  }

  sample.anno = sample.anno[!duplicated(Tumor_Sample_Barcode)] #sample.anno[Tumor_Sample_Barcode %in% maf.tsbs]
  anno.tsbs = sample.anno[,Tumor_Sample_Barcode]

  if(!length(maf.tsbs[!maf.tsbs %in% anno.tsbs]) == 0){
    if(chatty){
      cat('--Annotation missing for below samples in MAF:\n')
      for(temp in maf.tsbs[!maf.tsbs %in% anno.tsbs]){
        cat(paste0("  ", temp, "\n"))
      }
    }
  }
  sample.anno = sample.anno[Tumor_Sample_Barcode %in% maf.tsbs]

  return(list(variants.per.sample = tsb, variant.type.summary = vt.cast, variant.classification.summary = vc.cast,
              gene.summary = hs.cast, summary = summary, sample.anno = sample.anno))
}

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

flags = function(top = NULL){
  top100flags = c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", "MUC5B",
    "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2", "LAMA5", "AHNAK",
    "HMCN1", "USH2A", "DNAH11", "MACF1", "MUC17", "DNAH5", "GPR98",
    "FAT1", "PKD1", "MDN1", "RNF213", "RYR1", "DNAH2", "DNAH3", "DNAH8",
    "DNAH1", "DNAH9", "ABCA13", "APOB", "SRRM2", "CUBN", "SPTBN5",
    "PKHD1", "LRP2", "FBN3", "CDH23", "DNAH10", "FAT4", "RYR3", "PKHD1L1",
    "FAT2", "CSMD1", "PCNT", "COL6A3", "FRAS1", "FCGBP", "DNAH7",
    "RP1L1", "PCLO", "ZFHX3", "COL7A1", "LRP1B", "FAT3", "EPPK1",
    "VPS13C", "HRNR", "MKI67", "MYO15A", "STAB1", "ZAN", "UBR4",
    "VPS13B", "LAMA1", "XIRP2", "BSN", "KMT2C", "ALMS1", "CELSR1",
    "TG", "LAMA3", "DYNC2H1", "KMT2D", "BRCA2", "CMYA5", "SACS",
    "STAB2", "AKAP13", "UTRN", "VWF", "VPS13D", "ANK3", "FREM2",
    "PKD1L1", "LAMA2", "ABCA7", "LRP1", "ASPM", "MYOM2", "PDE4DIP",
    "TACC2", "MUC2", "TEP1", "HELZ2", "HERC2", "ABCA4")

  if(is.null(top)){
    top100flags
  }else{
    top100flags[1:top]
  }
}

loci2df = function(loci){
  chr = as.character(unlist(data.table::tstrsplit(x = loci, spli = ":", keep = 1)))
  start = unlist(data.table::tstrsplit(x = unlist(data.table::tstrsplit(x = loci, split = ":", keep = 2)), split = "-", keep = 1))
  start = as.numeric(as.character(gsub(pattern = ",", replacement = "", x = as.character(start))))
  end = unlist(data.table::tstrsplit(x = unlist(data.table::tstrsplit(x = loci, split = ":", keep = 2)), split = "-", keep = 2))
  end = as.numeric(as.character(gsub(pattern = ",", replacement = "", x = as.character(end))))
  data.table::data.table(Chromosome = chr, Start_Position = start, End_Position = end)
}
