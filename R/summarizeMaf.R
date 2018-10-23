summarizeMaf = function(maf, anno = NULL, chatty = TRUE){

  if('NCBI_Build' %in% colnames(maf)){
    NCBI_Build = unique(maf[!Variant_Type %in% 'CNV', NCBI_Build])
    NCBI_Build = NCBI_Build[!is.na(NCBI_Build)]

    if(chatty){
      if(length(NCBI_Build) > 1){
        message('NOTE: Mutiple reference builds found!')
        NCBI_Build = do.call(paste, c(as.list(NCBI_Build), sep=";"))
        message(NCBI_Build)
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
        message('Mutiple centers found.')
        print(Center)
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
  flags = c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", "MUC5B",
            "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2", "LAMA5", "AHNAK",
            "HMCN1", "USH2A", "DNAH11", "MACF1", "MUC17")

  #Variants per TSB
  tsb = maf[,.N, Tumor_Sample_Barcode]
  colnames(tsb)[2] = 'Variants'
  tsb = tsb[order(tsb$Variants, decreasing = TRUE),]

  #summarise and casting by 'Variant_Classification'
  vc = maf[,.N, .(Tumor_Sample_Barcode, Variant_Classification )]
  vc.cast = data.table::dcast(data = vc, formula = Tumor_Sample_Barcode ~ Variant_Classification, fill = 0, value.var = 'N')

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
  vt.cast = data.table::dcast(data = vt, formula = Tumor_Sample_Barcode ~ Variant_Type, value.var = 'N', fill = 0)

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

  if(chatty){
    print(summary)

    message("Gene Summary..")
    print(hs.cast)
  }

  #Check for flags.
  if(nrow(hs.cast) > 10){
    topten = as.character(hs.cast[1:10, Hugo_Symbol])
    topten = topten[topten %in% flags]
      if(chatty){
        if(length(topten) > 0){
          message('NOTE: Possible FLAGS among top ten genes:')
          print(topten)
      }
    }
  }


  if(chatty){
    message("Checking clinical data..")
  }

  if(is.null(anno)){
    if(chatty){
      message("NOTE: Missing clinical data! It is strongly recommended to provide clinical data associated with samples if available.")
    }
    sample.anno = data.table::data.table(Tumor_Sample_Barcode = maf.tsbs)
  }else if(is.data.frame(x = anno)){
      sample.anno  = data.table::setDT(anno)
      if(!'Tumor_Sample_Barcode' %in% colnames(sample.anno)){
        message(paste0('Available fields in provided annotations..'))
        print(colnames(sample.anno))
        stop(paste0('Tumor_Sample_Barcode column not found in provided clinical data. Rename column name containing sample names to Tumor_Sample_Barcode if necessary.'))
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
  sample.anno = as.data.frame(apply(sample.anno, 2, function(y) trimws(y))) #remove trailing whitespaces
  sample.anno[sample.anno == ""] = NA #Replace blanks with NA
  sample.anno = as.data.frame(apply(sample.anno, 2, function(y) gsub(pattern = " ", replacement = "_", x = y))) #replace spaces with _
  data.table::setDT(x = sample.anno)
  if(ncol(sample.anno) == 1){
    colnames(sample.anno)[1] = c("Tumor_Sample_Barcode")
  }

  sample.anno = sample.anno[!duplicated(Tumor_Sample_Barcode)] #sample.anno[Tumor_Sample_Barcode %in% maf.tsbs]
  anno.tsbs = sample.anno[,Tumor_Sample_Barcode]

  if(!length(maf.tsbs[!maf.tsbs %in% anno.tsbs]) == 0){
    if(chatty){
      message('Annotation missing for below samples in MAF')
      print(maf.tsbs[!maf.tsbs %in% anno.tsbs])
    }
  }
  sample.anno = sample.anno[Tumor_Sample_Barcode %in% maf.tsbs]

  return(list(variants.per.sample = tsb, variant.type.summary = vt.cast, variant.classification.summary = vc.cast,
              gene.summary = hs.cast, summary = summary, sample.anno = sample.anno))
}

get_vcColors = function(alpha = 1){
  col = c(RColorBrewer::brewer.pal(12,name = "Paired"), RColorBrewer::brewer.pal(11,name = "Spectral")[1:3],'black', 'violet', 'royalblue')
  col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  names(col) = names = c('Nonstop_Mutation','Frame_Shift_Del','IGR','Missense_Mutation','Silent','Nonsense_Mutation',
                         'RNA','Splice_Site','Intron','Frame_Shift_Ins','Nonstop_Mutation','In_Frame_Del','ITD','In_Frame_Ins',
                         'Translation_Start_Site',"Multi_Hit", 'Amp', 'Del')
  col
}

get_titvCol = function(alpha = 1){
  col = c('coral4', 'lightcyan4', 'cornflowerblue', 'lightsalmon1', 'forestgreen', 'deeppink3')
  col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  names(col) = c('C>T', 'C>G', 'C>A', 'T>A', 'T>C', 'T>G')
  col
}

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}
