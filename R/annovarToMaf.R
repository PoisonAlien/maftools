#' Converts annovar annotations into MAF.
#'
#' @description Converts variant annotations from Annovar into a basic MAF.
#' @details Annovar is one of the most widely used Variant Annotation tools in Genomics. Annovar output is generally in a tabular format with various annotation columns.
#' This function converts such annovar output files into MAF. This function requires that annovar was run with gene based annotation as a first operation, before including
#' any filter or region based annotations. Please be aware that this function performs no transcript prioritization.
#'
#' e.g,
#' table_annovar.pl example/ex1.avinput humandb/ -buildver hg19 -out myanno -remove -protocol (\code{refGene}),cytoBand,dbnsfp30a -operation (\code{g}),r,f -nastring NA
#'
#' This function mainly uses gene based annotations for processing, rest of the annotation columns from input file will be attached to the end of the resulting MAF.
#' @param annovar input annovar annotation file. Can be vector of multiple files.
#' @param Center Center field in MAF file will be filled with this value. Default NA.
#' @param refBuild NCBI_Build field in MAF file will be filled with this value. Default hg19.
#' @param tsbCol column name containing Tumor_Sample_Barcode or sample names in input file.
#' @param table reference table used for gene-based annotations. Can be 'ensGene' or 'refGene'. Default 'refGene'
#' @param basename If provided writes resulting MAF file to an output file.
#' @param sep field seperator for input file. Default tab seperated.
#' @param MAFobj If TRUE, returns results as an \code{\link{MAF}} object.
#' @param sampleAnno annotations associated with each sample/Tumor_Sample_Barcode in input annovar file. If provided it will be included in MAF object. Could be a text file or a data.frame. Ideally annotation would contain clinical data, survival information and other necessary features associated with samples. Default NULL.
#' @references Wang, K., Li, M. & Hakonarson, H. ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. Nucleic Acids Res 38, e164 (2010).
#' @return MAF table.
#' @examples
#' var.annovar <- system.file("extdata", "variants.hg19_multianno.txt", package = "maftools")
#' var.annovar.maf <- annovarToMaf(annovar = var.annovar, Center = 'CSI-NUS', refBuild = 'hg19',
#' tsbCol = 'Tumor_Sample_Barcode', table = 'ensGene')
#' @export

annovarToMaf = function(annovar, Center = NULL, refBuild = 'hg19', tsbCol = NULL, table = 'refGene', basename = NULL , sep = '\t', MAFobj = FALSE, sampleAnno = NULL){


  ann = lapply(annovar, data.table::fread, colClasses = 'character', sep = sep, stringsAsFactors = FALSE, fill = TRUE, header=TRUE)
  names(ann) = gsub(pattern = ".hg19_multianno.txt", replacement = "", x = basename(path = annovar))
  ann = data.table::rbindlist(l = ann, fill = TRUE, idcol = "sample_id")

  #Check to see if input file contains sample names
  if(is.null(tsbCol)){
    if(! 'Tumor_Sample_Barcode' %in% colnames(ann)){
      colnames(ann)[which(colnames(ann) == "sample_id")] = 'Tumor_Sample_Barcode'
      message("Tumor_Sample_Barcode field not found in input file. Using file names as Tumor_Sample_Barcode\nYou can manually specify column containing sample ids with argument `tsbCol`")
    }
  }else{
    colnames(ann)[which(colnames(ann) == tsbCol)] = 'Tumor_Sample_Barcode'
  }

  #Table options. See here: http://annovar.openbioinformatics.org/en/latest/user-guide/download/ (not considering UCSC known genes options for now)
  tabl.opts = c('refGene', 'ensGene')

  if(length(table) > 1){
    stop('table can only be either refGene or ensGene')
  }

  if(!table %in% tabl.opts){
    stop('table can only be either refGene or ensGene')
  }

  if(table == 'ensGene'){
    colnames(ann)[which(colnames(ann) == 'Func.ensGene')] = 'Func.refGene'
    colnames(ann)[which(colnames(ann) == 'Gene.ensGene')] = 'Gene.refGene'
    colnames(ann)[which(colnames(ann) == 'ExonicFunc.ensGene')] = 'ExonicFunc.refGene'
    colnames(ann)[which(colnames(ann) == 'AAChange.ensGene')] = 'AAChange.refGene'
    colnames(ann)[which(colnames(ann) == 'GeneDetail.ensGene')] = 'GeneDetail.refGene'
  }

    essential.col = c('Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene',
                      'ExonicFunc.refGene', 'AAChange.refGene')

#     #Change column names to standard names;
    for(i in 1:length(essential.col)){
      colId = suppressWarnings(grep(pattern = paste0('^', essential.col[i], '$'), x = colnames(ann), ignore.case = TRUE))
      if(length(colId) == 1){
        colnames(ann)[colId] = essential.col[i]
      }
     }

    if(length(essential.col[!essential.col %in% colnames(ann)]) > 0) {
      message('Available fields:')
      print(colnames(ann))
      message(paste0('Missing required field in input file: '))
      print(essential.col[!essential.col %in% colnames(ann)])
      stop()
    }


  #Center
  if(is.null(Center)){
    Center = NA
  }
  #Add unique ID for each variant
  ann$uid = paste('uid', 1:nrow(ann), sep='')

  #Mandatory fields
  ann.mand = c('Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'Tumor_Sample_Barcode', 'uid')

  #Rest of the optional fields (later they will be attached to the maf file)
  ann.opt = colnames(ann)[!colnames(ann) %in% ann.mand]
  ann.opt = c(ann.opt, 'uid')
  ann.opt = ann[,ann.opt,with = FALSE]

  ann = ann[,ann.mand,with = FALSE]

  ann$ExonicFunc.refGene = gsub(pattern = ' SNV', replacement = '', x = ann$ExonicFunc.refGene)

  funcSpl = strsplit(x = as.character(ann$ExonicFunc.refGene), split = ';', fixed = TRUE)
  funcSpl = sapply(funcSpl, function(l){l[length(l)]})
  ann$ExonicFunc.refGene = funcSpl

  funcRef = strsplit(x = as.character(ann$Func.refGene), split = ';', fixed = TRUE)
  funcRef = sapply(funcRef, function(l){l[length(l)]})
  ann$Func.refGene = funcRef

  #Change Variant Classification factors.
  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene == 'intronic', yes = 'Intron', no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene == 'intergenic', yes = 'IGR', no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene == 'downstream', yes = "3'Flank", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene == 'upstream', yes = "5'Flank", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene == 'splicing', yes = "Splice_Site", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene == 'UTR3', yes = "3'UTR", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene == 'UTR5', yes = "5'UTR", no = ann$ExonicFunc.refGene)

  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene %in% c('ncRNA_exonic', 'ncRNA_intronic', 'ncRNA_UTR3', 'ncRNA_UTR5', 'ncRNA'),
                                  yes = 'RNA', no = ann$ExonicFunc.refGene)

  ann.lvls = c('synonymous', 'nonsynonymous', 'stopgain', 'stoploss', 'frameshift insertion', 'frameshift deletion', 'nonframeshift insertion',
               'nonframeshift deletion', 'Intron', 'IGR', 'Splice_Site', "3'UTR", "3'Flank", "5'UTR", "5'Flank", "unknown", "UNKNOWN", 'RNA')
               # "intronic", "intergenic", "downstream", "upstream", "splicing", "UTR3", "UTR5",
               # 'ncRNA_exonic', 'ncRNA_intronic', 'ncRNA_UTR3', 'ncRNA_UTR5', 'ncRNA')
  ann.lbls = c('Silent', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Frame_Shift_Ins', 'Frame_Shift_Del', 'In_Frame_Ins',
               'In_Frame_Del', 'Intron', 'IGR', 'Splice_Site', "3'UTR", "3'Flank", "5'UTR", "5'Flank", "UNKNOWN", "UNKNOWN", 'RNA')
               # "Intron", "'IGR", "3'Flank", "5'Flank", "Splice_Site", "3'UTR", "5'UTR",
               # 'RNA', 'RNA', 'RNA', 'RNA', 'RNA')
  names(ann.lbls) = ann.lvls

  ann$ExonicFunc.refGene = as.character(ann.lbls[as.character(ann$ExonicFunc.refGene)])

  #ann$ExonicFunc.refGene = suppressWarnings(as.character( factor(x = ann$ExonicFunc.refGene, levels = ann.lvls, labels = ann.lbls)))

  #Chnage the way indels are representaed.
  ann.del = ann[ann$Alt %in% "-"]
  ann = ann[!ann$Alt %in% "-"]

  if(nrow(ann.del) > 0){
    ann.del$var.type = 'DEL'
  }

  ann.ins = ann[ann$Ref %in% "-"]
  ann = ann[!ann$Ref %in% "-"]

  if(nrow(ann.ins) > 0){
    ann.ins$var.type = 'INS'
  }

  if(nrow(ann) > 0){
    ann$var.type = 'SNP'
  }

  ann = rbind(ann, ann.del, ann.ins, fill = TRUE)

  ann.splice = ann[ExonicFunc.refGene == 'Splice_Site']
  if(nrow(ann.splice) > 0){
    ann = ann[ExonicFunc.refGene != 'Splice_Site']
    #ann.splice$AAChange.refGene = gsub(x = sapply(strsplit(x = as.character(ann.splice$Gene.refGene), split = '(', fixed = T), '[[',2), pattern = ')$', replacement = '')
    ann.splice$Gene.refGene = sapply(strsplit(x = as.character(ann.splice$Gene.refGene), split = '(', fixed = TRUE), '[[',1)
    ann = rbind(ann, ann.splice, fill = TRUE)
  }

  #protein and tx changes
  #NOTE: for now last transcript is considered from the total annotaion.
  xaa = strsplit(as.character(ann$AAChange.refGene),split = ':',fixed = TRUE)
  proteinChange = sapply(xaa, function(l){l[length(l)]})
  tx = unlist(sapply(xaa, function(l){l[length(l)-3]}))
  txChange = unlist(sapply(xaa, function(l){if(length(l) > 1){l[length(l)-1]} else{NA}}))

  #Make final maf table
  ann.maf = data.table::data.table(Hugo_Symbol = ann$Gene.refGene, Entrez_Gene_Id = NA, Center = Center, NCBI_Build = refBuild, Chromosome = ann$Chr, Start_Position = ann$Start, End_Position = ann$End, Strand = '+',
                       Variant_Classification = ann$ExonicFunc.refGene, Variant_Type = ann$var.type, Reference_Allele = ann$Ref, Tumor_Seq_Allele1 = ann$Ref, Tumor_Seq_Allele2 = ann$Alt,
                       dbSNP_RS = NA, Tumor_Sample_Barcode = ann$Tumor_Sample_Barcode, Mutation_Status = 'Somatic',
                       AAChange = proteinChange, Transcript_Id = tx, TxChange = txChange, uid = ann$uid)

  ann.maf = merge(ann.maf, ann.opt, by = 'uid')
  ann.maf = ann.maf[,uid := NULL] #Remove unique ids.

  #Annovar ensGene doesn't provide HGNC gene symbols as Hugo_Symbol. We will change them manually.
  if(table == 'ensGene'){
    ens = system.file('extdata', 'ensGenes.txt.gz', package = 'maftools')
    message('Converting Ensemble Gene IDs into HGNC gene symbols.')
    if(Sys.info()[['sysname']] == 'Windows'){
      ens.gz = gzfile(description = ens, open = 'r')
      ens <- suppressWarnings( data.table(read.csv( file = ens.gz, header = TRUE, sep = '\t', stringsAsFactors = FALSE)) )
      close(ens.gz)
    } else{
      ens = data.table::fread(cmd = paste('zcat <', ens), sep = '\t', stringsAsFactors = FALSE)
    }

    ann.maf = merge(ann.maf, ens, by.x = 'Hugo_Symbol', by.y = 'ens_id', all.x = TRUE)
    ann.maf[,ens_id := Hugo_Symbol] #Backup original ids
    ann.maf[,Hugo_Symbol := hgnc_symbol] #Add GHNC gene names
    ann.maf[,Entrez_Gene_Id := Entrez] #Add entrez identifiers.
    message('Done! Original ensemble gene IDs are preserved under field name ens_id')
  }

  if(!is.null(basename)){
    write.table(x = ann.maf, file = paste(basename, 'maf', sep='.'), sep='\t', quote = FALSE, row.names = FALSE)
  }

  if(MAFobj){
    #Convert to factors.
    ann.maf$Tumor_Sample_Barcode = as.factor(x = as.character(ann.maf$Tumor_Sample_Barcode))
    ann.maf$Variant_Classification = as.factor(x = as.character(ann.maf$Variant_Classification))
    ann.maf$Variant_Type = as.factor(x = as.character(ann.maf$Variant_Type))
    ann.maf = validateMaf(maf = ann.maf, isTCGA = FALSE, rdup = TRUE)
    ann.maf.summary = summarizeMaf(maf = ann.maf, anno = sampleAnno, chatty = FALSE)
    if(length(unique(ann.maf[,Tumor_Sample_Barcode])) < 2){
      message('Too few samples to create MAF object. Returning MAF table.')
      return(ann.maf)
    }else{
      silent = c("3'UTR", "5'UTR", "3'Flank", "Targeted_Region", "Silent", "Intron",
                 "RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA")
      ann.maf.silent = ann.maf[Variant_Classification %in% silent]
      ann.maf = ann.maf[!Variant_Classification %in% silent]

      m = MAF(data = ann.maf, variants.per.sample = ann.maf.summary$variants.per.sample, variant.type.summary = ann.maf.summary$variant.type.summary,
              variant.classification.summary = ann.maf.summary$variant.classification.summary, gene.summary = ann.maf.summary$gene.summary,
              summary = ann.maf.summary$summary, maf.silent = ann.maf.silent, clinical.data = ann.maf.summary$sample.anno)

      return(m)
    }

  }else{
    # For reasons returned data.table doesn't print in first attempt.
    # Discussion here: https://github.com/Rdatatable/data.table/issues/939
    # As a solution try to print before returning.
    print(ann.maf)
    return(ann.maf)
  }

}
