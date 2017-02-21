validateMaf = function(maf, rdup = TRUE, isTCGA = isTCGA, chatty = TRUE){

  #necessary fields.
  required.fields = c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
                      'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode')

  #Change column names to standard names; i.e, camel case
  for(i in 1:length(required.fields)){
    colId = suppressWarnings(grep(pattern = required.fields[i], x = colnames(maf), ignore.case = TRUE))
    if(length(colId) > 0){
      colnames(maf)[colId] = required.fields[i]
    }
  }

  missing.fileds = required.fields[!required.fields %in% colnames(maf)] #check if any of them are missing

  if(length(missing.fileds) > 0){
    missing.fileds = paste(missing.fileds[1], sep = ',', collapse = ', ')
    stop(paste('missing required fields from MAF:', missing.fileds)) #stop if any of required.fields are missing
  }

  #convert "-" to "." in "Tumor_Sample_Barcode" to avoid complexity in naming
  maf$Tumor_Sample_Barcode = gsub(pattern = '-', replacement = '.', x = as.character(maf$Tumor_Sample_Barcode))

  if(rdup){
    maf = maf[, variantId := paste(Chromosome, Start_Position, Tumor_Sample_Barcode, sep = ':')]
    if(nrow(maf[duplicated(variantId)]) > 0){
      if(chatty){
        message("NOTE: Removed ",  nrow(maf[duplicated(variantId)]) ," duplicated variants")
      }
      maf = maf[!duplicated(variantId)]
    }
    maf[,variantId := NULL]
  }

  if(nrow(maf[Hugo_Symbol %in% ""]) > 0){
    if(chatty){
      message('NOTE: Found ', nrow(maf[Hugo_Symbol %in% ""]), ' variants with no Gene Symbols.')
      print(maf[Hugo_Symbol %in% "", required.fields, with = FALSE])
      message("Annotating them as 'UnknownGene' for convenience")
    }
    maf$Hugo_Symbol = ifelse(test = maf$Hugo_Symbol == "", yes = 'UnknownGene', no = maf$Hugo_Symbol)
  }

  if(nrow(maf[is.na(Hugo_Symbol)]) > 0){
    if(chatty){
      message('NOTE: Found ', nrow(maf[is.na(Hugo_Symbol) > 0]), ' variants with no Gene Symbols.')
      print(maf[is.na(Hugo_Symbol), required.fields, with =FALSE])
      message("Annotating them as 'UnknownGene' for convenience")
    }
    maf$Hugo_Symbol = ifelse(test = is.na(maf$Hugo_Symbol), yes = 'UnknownGene', no = maf$Hugo_Symbol)
  }

  if(isTCGA){
    maf$Tumor_Sample_Barcode = substr(x = maf$Tumor_Sample_Barcode, start = 1, stop = 12)
  }

  #Variant Classification with Low/Modifier variant consequences. http://asia.ensembl.org/Help/Glossary?id=535
  silent = c("3'UTR", "5'UTR", "3'Flank", "Targeted_Region", "Silent", "Intron",
             "RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA")
  #Variant Classification with High/Moderate variant consequences. http://asia.ensembl.org/Help/Glossary?id=535
  vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                   "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                   "In_Frame_Ins", "Missense_Mutation")
  vt = c('SNP', 'DNP', 'TNP', 'ONP', 'INS', 'DEL')

  maf.vcs = unique(as.character(maf[,Variant_Classification]))
  maf.vts = unique(as.character(maf[,Variant_Type]))

  if(length(maf.vcs[!maf.vcs %in% c(silent, vc.nonSilent)] > 0)){
    message("NOTE: Non MAF specific values in Variant_Classification column:")
    print(maf.vcs[!maf.vcs %in% c(silent, vc.nonSilent)])
  }

  if(length(maf.vts[!maf.vts %in% vt] > 0)){
    message("NOTE: Non MAF specific values in Variant_Type column:")
    print(maf.vts[!maf.vts %in% vt])
  }

  return(maf)
}
