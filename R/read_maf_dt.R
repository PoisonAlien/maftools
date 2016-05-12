#' Read MAF files.
#' @description Takes tab delimited MAF (can be plain text or gz compressed) file as an input and summarizes it in various ways. Also creates oncomatrix - helpful for visualization.
#'
#' @param maf tab delimited MAF file. File can also be gz compressed.
#' @param removeSilent logical. Whether to discard silent (variants with Low/Modifier consequences) mutations ("3'UTR", "5'UTR", "3'Flank", "Targeted_Region", "Silent", "Intron","RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA"). Default is TRUE.
#' @param useAll logical. Whether to use all variants irrespective of values in Mutation_Status. Defaults to False. Only uses with values Somatic.
#' @return An object of class MAF.
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf, removeSilent = TRUE, useAll = FALSE)
#' @import data.table
#' @import ggplot2
#' @seealso \code{\link{plotmafSummary}} \code{\link{write.mafSummary}}
#' @export


read.maf = function(maf, removeSilent = TRUE, useAll = FALSE){

  message('reading maf..')

  if(as.logical(length(grep(pattern = 'gz$', x = maf, fixed = FALSE)))){
    #If system is Linux use fread, else use gz connection to read gz file.
    if(Sys.info()[['sysname']] == 'Windows'){
      maf.gz = gzfile(description = maf, open = 'r')
      suppressWarnings(maf <- data.table(read.csv(file = maf.gz, header = TRUE, sep = '\t', stringsAsFactors = FALSE)))
      close(maf.gz)
    } else{
      maf = suppressWarnings(fread(input = paste('zcat <', maf), sep = '\t', stringsAsFactors = FALSE, verbose = FALSE, data.table = TRUE, showProgress = TRUE, header = TRUE))
    }
  } else{
    suppressWarnings(maf <- fread(input = maf, sep = "\t", stringsAsFactors = FALSE, verbose = FALSE, data.table = TRUE, showProgress = TRUE, header = TRUE))
  }

  #validate MAF file
  maf = validateMaf(maf = maf)

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

  #Variant Classification with Low/Modifier variant consequences. http://asia.ensembl.org/Help/Glossary?id=535
  silent = c("3'UTR", "5'UTR", "3'Flank", "Targeted_Region", "Silent", "Intron",
             "RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA")
  #Variant Classification with High/Moderate variant consequences. http://asia.ensembl.org/Help/Glossary?id=535
  vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                   "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                   "In_Frame_Ins", "Missense_Mutation")

  maf.silent = maf[Variant_Classification %in% silent]

  if(removeSilent){

    if(nrow(maf.silent) > 0){
      maf.silent.vc = maf.silent[,.N, .(Tumor_Sample_Barcode, Variant_Classification)]
      maf.silent.vc.cast = data.table::dcast(data = maf.silent.vc, formula = Tumor_Sample_Barcode ~ Variant_Classification, fill = 0, value.var = 'N') #why dcast is not returning it as data.table ?
      summary.silent = data.table(ID = c('Samples',colnames(maf.silent.vc.cast)[2:ncol(maf.silent.vc.cast)]),
                                  N = c(nrow(maf.silent.vc.cast), colSums(maf.silent.vc.cast[,2:ncol(maf.silent.vc.cast), with = FALSE])))

      maf = maf[!Variant_Classification %in% silent] #Remove silent variants from main table
      message(paste('Excluding',nrow(maf.silent), 'silent variants.'))
      print(summary.silent)
    } else{
      message(message(paste('Excluding',nrow(maf.silent), 'silent variants.')))
    }
  }else{
    message('Silent variants are being kept!')
  }

  #convert to factors
  maf$Variant_Type = as.factor(as.character(maf$Variant_Type))
  maf$Variant_Classification = as.factor(as.character(maf$Variant_Classification))
  maf$Tumor_Sample_Barcode = as.factor(as.character(maf$Tumor_Sample_Barcode))

  message('Summarizing..')
  mafSummary = summarizeMaf(maf = maf)
  print(mafSummary$summary)

  message("Frequently mutated genes..")
  print(mafSummary$gene.summary)

  #Create oncomatrix
  oncomat = createOncoMatrix(maf)

  if(is.null(oncomat)){
    m = MAF(data = maf, variants.per.sample = mafSummary$variants.per.sample, variant.type.summary = mafSummary$variant.type.summary,
            variant.classification.summary = mafSummary$variant.classification.summary,gene.summary = mafSummary$gene.summary,
            oncoMatrix = NULL, numericMatrix = NULL, summary = mafSummary$summary,
            classCode = NULL, maf.silent = maf.silent)
  }else{
    m = MAF(data = maf, variants.per.sample = mafSummary$variants.per.sample, variant.type.summary = mafSummary$variant.type.summary,
            variant.classification.summary = mafSummary$variant.classification.summary,gene.summary = mafSummary$gene.summary,
            oncoMatrix = oncomat$oncomat, numericMatrix = oncomat$nummat, summary = mafSummary$summary,
            classCode = oncomat$vc, maf.silent = maf.silent)
  }

  message('Done !')
  return(m)
}
