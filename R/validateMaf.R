validateMaf = function(maf, isTCGA = isTCGA){

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

  if(isTCGA){
    maf$Tumor_Sample_Barcode = substr(x = maf$Tumor_Sample_Barcode, start = 1, stop = 12)
  }

  return(maf)
}


