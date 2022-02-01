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
#' @param ens2hugo If `table` is `ensGene`, setting this argument to `TRUE` converts all ensemble IDs to hugo symbols.
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

annovarToMaf = function(annovar, Center = NULL, refBuild = 'hg19', tsbCol = NULL, table = 'refGene', ens2hugo = TRUE, basename = NULL , sep = '\t', MAFobj = FALSE, sampleAnno = NULL){

  start_time = proc.time()
  cat("-Reading annovar files\n")
  ann = lapply(annovar, data.table::fread, colClasses = 'character', sep = sep, stringsAsFactors = FALSE,
               fill = TRUE, header=TRUE, skip = "Chr")
  names(ann) = unlist(data.table::tstrsplit(x = basename(annovar), split = "\\.", keep = 1))
  ann = data.table::rbindlist(l = ann, fill = TRUE, idcol = "sample_id", use.names = TRUE)

  #Check to see if input file contains sample names
  if(is.null(tsbCol)){
    if(! 'Tumor_Sample_Barcode' %in% colnames(ann)){
      colnames(ann)[which(colnames(ann) == "sample_id")] = 'Tumor_Sample_Barcode'
      cat("--Tumor_Sample_Barcode column not found. Creating sample IDs from filenames\n")
    }
  }else{
    colnames(ann)[which(colnames(ann) == tsbCol)] = 'Tumor_Sample_Barcode'
  }

  #Table options. See here: http://annovar.openbioinformatics.org/en/latest/user-guide/download/ (not considering UCSC known genes options for now)
  table = match.arg(arg = table, choices = c('refGene', 'ensGene'))

  if(table == 'ensGene'){
    colnames(ann)[which(colnames(ann) == 'Func.ensGene')] = 'Func.refGene'
    colnames(ann)[which(colnames(ann) == 'Gene.ensGene')] = 'Gene.refGene'
    colnames(ann)[which(colnames(ann) == 'ExonicFunc.ensGene')] = 'ExonicFunc.refGene'
    colnames(ann)[which(colnames(ann) == 'AAChange.ensGene')] = 'AAChange.refGene'
    colnames(ann)[which(colnames(ann) == 'GeneDetail.ensGene')] = 'GeneDetail.refGene'
  }

    essential.col = c('Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene',
                      'ExonicFunc.refGene', 'AAChange.refGene')

    #Change column names to standard names;
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

    #If multiple genes are assigned, used the first entry (which correpsonds to closest gene)
    ann[, Hugo_Symbol := unlist(data.table::tstrsplit(Gene.refGene, split = ";", keep = 1))]

    #Annovar to MAF mappings (http://annovar.openbioinformatics.org/en/latest/user-guide/gene/)
    annovar_values = c(
      exonic = "RNA",
      splicing = "Splice_Site",
      UTR5 = "5'UTR",
      UTR3 = "3'UTR",
      intronic = "Intron",
      upstream = "5'Flank",
      downstream = "3'Flank",
      intergenic = "IGR",
      `frameshift insertion` = "Frame_Shift_Ins",
      `frameshift deletion` = "Frame_Shift_Del",
      `frameshift block substitution` = "Frameshift_INDEL",
      `frameshift substitution` = "Frameshift_INDEL",
      stopgain = "Nonsense_Mutation",
      stoploss = "Nonstop_Mutation",
      startloss = "Translation_Start_Site",
      startgain = "Unknown", #Can not properly map MAF classification. Setting it to Unknown
      `nonframeshift insertion` = "In_Frame_Ins",
      `nonframeshift deletion` = "In_Frame_Del",
      `nonframeshift block substitution` = "Inframe_INDEL",
      `nonframeshift substitution` = "Inframe_INDEL",
      `nonsynonymous SNV` = "Missense_Mutation",
      `synonymous SNV` = "Silent",
      unknown = "Unknown",
      ncRNA_exonic = "RNA",
      ncRNA_intronic = "RNA",
      ncRNA_UTR3 = "RNA",
      ncRNA_UTR5 = "RNA",
      ncRNA = "RNA",
      ncRNA_splicing = "RNA"
    )

    ann_exonic = ann[Func.refGene %in% 'exonic']
    ann_res = ann[!Func.refGene %in% 'exonic']
    if(nrow(ann_exonic) ==0 & nrow(ann_res) == 0){
      stop("No suitable exonic or intronic variants found!")
    }

    #Exonic variants
    if(nrow(ann_exonic) > 0){
      cat("-Processing Exonic variants\n")

      #Choose first entry for mixed annotations (e.g; exonic;splicing)
      ann_exonic[, Func.refGene := data.table::tstrsplit(x = as.character(ann_exonic$Func.refGene), split = ";", keep = 1)]
      cat("--Adding Variant_Classification\n")
      ann_exonic[,Variant_Classification := annovar_values[ExonicFunc.refGene]]

      #Parse aa-change for exonic variants
      cat("--Parsing aa-change\n")
      aa_change = unlist(data.table::tstrsplit(x = as.character(ann_exonic$AAChange.refGene),split = ',', fixed = TRUE, keep = 1))
      aa_tbl = lapply(aa_change, function(x){
        x = unlist(strsplit(x = x, split = ":", fixed = TRUE))

        if(length(x) == 5){
          #column contains the gene name, the transcript identifier, exon, tx-change, aa-change
          #If these entries are missing, fill them with NAs
          tx = x[2]
          exon = x[3]
          txChange = x[4]
          aaChange = x[5]
        }else{
          tx = NA
          exon = NA
          txChange = NA
          aaChange = NA
        }

        data.table::data.table(tx, exon, txChange, aaChange)
      })
      aa_tbl = data.table::rbindlist(l = aa_tbl)

      if(length(aa_change) != nrow(ann_exonic)){
        stop("Something went wrong parsing aa-change")
      }
      ann_exonic = cbind(ann_exonic, aa_tbl)
    }

    #Non-exonic variants
    if(nrow(ann_res) > 0){
      cat("-Processing Non-exonic variants\n")
      ann_res[, Func.refGene := data.table::tstrsplit(x = as.character(ann_res$Func.refGene), split = ";", keep = 1)]
      cat("--Adding Variant_Classification\n")
      ann_res[,Variant_Classification := annovar_values[Func.refGene]]
      #Merge exonic and non-exonic variants
      ann = data.table::rbindlist(l = list(ann_exonic, ann_res), use.names = TRUE, fill = TRUE)
    }else{
      ann = ann_exonic
    }

    #Add Variant-type annotations based on difference between ref and alt alleles
    cat("-Adding Variant_Type\n")
    ann[,ref_alt_diff := nchar(Ref) - nchar(Alt)]
    ann[, Variant_Type := ifelse(
      ref_alt_diff == 0 ,
      yes = ifelse(test = Ref != "-" & Alt != "-", yes = "SNP",
      no = ifelse(ref_alt_diff < 0 , yes = "INS", no = "DEL")), no = NA)
    ]

    #Check for MNPs (they are neither INDELS nor SNPs)
    ann$Variant_Type = ifelse(
        test = ann$Variant_Type == "SNP",
        yes = ifelse(
          test = nchar(ann$Ref) > 1,
          yes = "MNP",
          no = "SNP"
        ),
        no = ann$Variant_Type
      )

    #Annotate MNPs as Unknown VC
    ann_mnps = ann[Variant_Type %in% "MNP"]
    if(nrow(ann_mnps) > 0){
      ann = ann[!Variant_Type %in% "MNP"]
      ann_mnps[,Variant_Classification := "Unknown"]
      ann = rbind(ann, ann_mnps)
      rm(ann_mnps)
    }

    #Frameshift_INDEL, Inframe_INDEL are annotated by annovar without INS or DEL status
    #Parse this based on difference between ref and alt alleles
    ann_indel = ann[Variant_Classification %in% c("Frameshift_INDEL", "Inframe_INDEL")]
    if(nrow(ann_indel) > 0){
      cat("-Fixing ambiguous INDEL annotations\n")
      ann = ann[!Variant_Classification %in% c("Frameshift_INDEL", "Inframe_INDEL")]
      vc_fixed = lapply(1:nrow(ann_indel), function(i) {
        x = ann_indel[i, Variant_Classification]
        if (x == "Frameshift_INDEL") {
          if (ann_indel[i, ref_alt_diff] > 0) {
            return("Frame_Shift_Del")
          } else{
            return("Frame_Shift_Ins")
          }
        } else if (x == "Inframe_INDEL") {
          if (ann_indel[i, ref_alt_diff] > 0) {
            return("In_Frame_Del")
          } else{
            return("In_Frame_Ins")
          }
        } else{
          return(x)
        }
      })
      ann_indel[,Variant_Classification := unlist(vc_fixed)]
      ann = rbind(ann, ann_indel)
    }

  #Annovar ensGene doesn't provide HGNC gene symbols as Hugo_Symbol. We will change them manually.
  if(table == 'ensGene'){
    if(ens2hugo){
      ens = system.file('extdata', 'ensGenes.txt.gz', package = 'maftools')
      cat('-Converting Ensemble Gene IDs into HGNC gene symbols\n')
      ens = data.table::fread(file =  ens, sep = '\t', stringsAsFactors = FALSE)

      ann = merge(ann, ens, by.x = 'Hugo_Symbol', by.y = 'ens_id', all.x = TRUE)
      ann[,ens_id := Hugo_Symbol] #Backup original ids
      ann[,Hugo_Symbol := hgnc_symbol] #Add GHNC gene names
      ann[,Entrez_Gene_Id := Entrez] #Add entrez identifiers.
      cat('--Done. Original ensemble gene IDs are preserved under field name ens_id\n')
    }
  }
    ann[, ref_alt_diff := NULL]

  #Re-roganize columns
  colnames(ann)[match(c("Chr", "Start", "End", "Ref", "Alt"), colnames(ann))] = c("Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2")
  ord1 = colnames(x = ann)[colnames(x = ann) %in% c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode", "tx", "exon", "txChange", "aaChange")]
  ord2 = colnames(x = ann)[!colnames(x = ann) %in% c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode", "tx", "exon", "txChange", "aaChange")]
  ann = ann[,c(ord1, ord2), with = FALSE]


  if(!is.null(basename)){
    data.table::fwrite(x = ann, file = paste(basename, 'maf', sep='.'), sep='\t')
  }

  cat("Finished in",data.table::timetaken(start_time),"\n")

  if(MAFobj){
    m = read.maf(maf = ann, clinicalData = sampleAnno, verbose = FALSE)
    return(m)
  }else{
    return(ann)
  }
}
