#' Converts ICGC Simple Somatic Mutation format file to MAF
#'
#' @description Converts ICGC Simple Somatic Mutation format file to Mutation Annotation Format. Basic fields are converted as per MAF specififcations, rest of the fields are retained as in the input file.
#' Ensemble gene IDs are converted to HGNC Symbols. Note that by default Simple Somatic Mutation format contains all affected transcripts of a variant resuting in
#' multiple entries of the same variant in same sample. It is hard to choose a single affected transcript based on annotations alone and by default this program removes repeated variants as duplicated entries.
#' If you wish to keep all of them, set removeDuplicatedVariants to FALSE.
#' @details ICGC Simple Somatic Mutattion format specififcation can be found here: http://docs.icgc.org/submission/guide/icgc-simple-somatic-mutation-format/
#'
#' @param icgc Input data in ICGC Simple Somatic Mutation format. Can be gz compressed.
#' @param basename If given writes to output file with basename.
#' @param MAFobj If TRUE returns results as an \code{\link{MAF}} object.
#' @param removeDuplicatedVariants removes repeated variants in a particuar sample, mapped to multiple transcripts of same Gene. See Description. Default TRUE.
#' @param addHugoSymbol If TRUE replaces ensemble gene IDs with Hugo_Symbols. Default FALSE.
#' @return tab delimited MAF file.
#' @examples
#' esca.icgc <- system.file("extdata", "simple_somatic_mutation.open.ESCA-CN.sample.tsv.gz", package = "maftools")
#' esca.maf <- icgcSimpleMutationToMAF(icgc = esca.icgc)
#' @export

icgcSimpleMutationToMAF = function(icgc, basename = NA, MAFobj = FALSE, removeDuplicatedVariants = TRUE, addHugoSymbol = FALSE){

  if(as.logical(length(grep(pattern = 'gz$', x = icgc, fixed = FALSE)))){

    if(Sys.info()[['sysname']] == 'Windows'){
      icgc.gz = gzfile(description = icgc, open = 'r')
      icgc <- suppressWarnings( data.table(read.csv( file = icgc.gz, header = TRUE, sep = '\t', stringsAsFactors = FALSE)) )
      close(icgc.gz)
    }else{
      icgc = suppressWarnings(data.table::fread(input = paste('zcat <', icgc), sep = '\t', stringsAsFactors = FALSE, verbose = FALSE, data.table = TRUE, showProgress = TRUE, header = TRUE))
    }
  }else{
    icgc = data.table::fread(input = icgc, sep = '\t', stringsAsFactors = FALSE, header = TRUE, showProgress = TRUE, data.table = TRUE, verbose = FALSE)
  }

  icgc.indels = icgc[consequence_type %in% 'frameshift_variant']
  icgc = icgc[!consequence_type %in% 'frameshift_variant']


  #Variant consequences
  vep.vc = c("splice_acceptor_variant", "splice_donor_variant", "transcript_ablation",
    "exon_loss_variant", "stop_gained", "stop_lost", "initiator_codon_variant",
    "start_lost", "missense_variant", "coding_sequence_variant",
    "conservative_missense_variant", "rare_amino_acid_variant", "transcript_amplification",
    "intron_variant", "INTRAGENIC", "intragenic_variant", "splice_region_variant",
    "incomplete_terminal_codon_variant", "synonymous_variant", "stop_retained_variant",
    "NMD_transcript_variant", "mature_miRNA_variant", "exon_variant",
    "non_coding_exon_variant", "non_coding_transcript_exon_variant",
    "non_coding_transcript_variant", "nc_transcript_variant", "5_prime_UTR_variant",
    "5_prime_UTR_premature_start_codon_gain_variant", "3_prime_UTR_variant",
    "TF_binding_site_variant", "regulatory_region_variant", "regulatory_region",
    "intergenic_variant", "intergenic_region", "upstream_gene_variant",
    "downstream_gene_variant", "disruptive_inframe_deletion", "inframe_deletion", "inframe_insertion", "disruptive_inframe_insertion")

  #Corresponding MAF Variant Classifications
  maf.vc = c("Splice_Site", "Splice_Site", "Splice_Site", "Splice_Site",
    "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site",
    "Translation_Start_Site", "Missense_Mutation", "Missense_Mutation",
    "Missense_Mutation", "Missense_Mutation", "Intron", "Intron",
    "Intron", "Intron", "Splice_Region", "Silent", "Silent", "Silent",
    "Silent", "RNA", "RNA", "RNA", "RNA", "RNA", "RNA", "5'UTR",
    "5'UTR", "3'UTR", "IGR", "IGR", "IGR", "IGR", "IGR", "5'Flank",
    "3'Flank", "In_Frame_Del", "In_Frame_Del", "In_Frame_Ins", "In_Frame_Ins")

  names(maf.vc) = vep.vc

  #Add MAF Variant_Classification levels correspoding to VEP consequences
  icgc$Variant_Classification = maf.vc[as.character(icgc$consequence_type)]
  #icgc$Variant_Classification = suppressWarnings(as.character(factor(x = icgc$consequence_type, levels = vep.vc, labels = maf.vc)))

  #Add MAF Variant_Classification levels of INDELS correspoding to VEP consequences
  icgc.indels$Variant_Classification = ifelse(test = icgc.indels$consequence_type %in% "frameshift_variant" & icgc.indels$mutation_type %in% "insertion of <=200bp", yes = "Frame_Shift_Ins", no = "Frame_Shift_Del")

  #bind everything
  icgc = rbind(icgc, icgc.indels)

  #Add MAF Variant_Type to correspoding VEP mutation_type
  icgc$Variant_Type = suppressWarnings(as.character(factor(x = icgc$mutation_type, levels = c("deletion of <=200bp", "insertion of <=200bp", "single base substitution"),
                             labels = c("DEL", "INS", "SNP"))))

  icgc.fields = colnames(icgc)

  required.fields = c("gene_affected", "assembly_version", "chromosome", "chromosome_start", "chromosome_end", "chromosome_strand", "Variant_Classification", "Variant_Type",
    "reference_genome_allele", "mutated_from_allele", "mutated_to_allele", "icgc_sample_id", "verification_status", "biological_validation_status", "sequencing_strategy",
    "verification_platform")

  icgc.rest.fields = icgc.fields[!icgc.fields %in% required.fields]

  icgc.rest = icgc[,icgc.rest.fields, with = FALSE]

  #Make MAF table
  icgc.maf = data.table::data.table(Hugo_Symbol = icgc[,gene_affected], Entrez_Gene_Id = NA, Center = NA, NCBI_Build = icgc[,assembly_version], Chromosome = icgc[,chromosome],
                         Start_Position = icgc[,chromosome_start], End_Position = icgc[,chromosome_end], Strand = '+',
                         Variant_Classification = icgc[,Variant_Classification], Variant_Type = icgc[,Variant_Type], Reference_Allele = icgc[,reference_genome_allele],
                         Tumor_Seq_Allele1 = icgc[,mutated_from_allele], Tumor_Seq_Allele2 = icgc[,mutated_to_allele], dbSNP_RS = NA, dbSNP_Val_Status  = NA,
                         Tumor_Sample_Barcode = icgc[,icgc_sample_id], Verification_Status = icgc[,verification_status], Mutation_Status = "Somatic",
                         Sequence_Source = icgc[,sequencing_strategy], Validation_Method = icgc[,verification_platform])

  icgc.maf$Verification_Status = ifelse(test = icgc.maf$Verification_Status %in% "not tested", yes = 'Unknown', no = 'Verified')


  #Attach rest of fields
  icgc.maf = cbind(icgc.maf, icgc.rest)

  if(addHugoSymbol){

    #Change ensemble gene IDs into Hugo_Symbol
    ens = system.file('extdata', 'ensGenes.txt.gz', package = 'maftools')
    message('Converting Ensemble Gene IDs into HGNC gene symbols.')

    if(Sys.info()[['sysname']] == 'Windows'){
      ens.gz = gzfile(description = ens, open = 'r')
      ens <- suppressWarnings( data.table(read.csv( file = ens.gz, header = TRUE, sep = '\t', stringsAsFactors = FALSE)) )
      close(ens.gz)
    } else{
      ens = data.table::fread(input = paste('zcat <', ens), sep = '\t', stringsAsFactors = FALSE)
    }

    icgc.maf[,ens_id := Hugo_Symbol]
    icgc.maf = merge(icgc.maf, ens, by.x = 'Hugo_Symbol', by.y = 'ens_id', all.x = TRUE)

    icgc.maf[,Hugo_Symbol := hgnc_symbol]
    icgc.maf[,Entrez_Gene_Id := Entrez]
    icgc.maf[,hgnc_symbol := NULL]
    icgc.maf[,Entrez := NULL]

    message('Done! Original ensemble gene IDs are preserved under field name ens_id')
  }


  #Order according to Hugo_Symbol
  icgc.maf = icgc.maf[order(Hugo_Symbol)]

  #Validate
  icgc.maf = validateMaf(maf = icgc.maf, isTCGA = FALSE, rdup = removeDuplicatedVariants)

  if(!is.na(basename)){
    write.table(x = icgc.maf, file = paste(basename, 'maf', sep='.'), quote = FALSE, row.names = FALSE, sep = '\t')
  }

  if(MAFobj){

    if(length(unique(icgc.maf[,Tumor_Sample_Barcode])) < 2){
      message('Too few samples to create MAF object. Returning MAF table.')
      return(icgc.maf)
    }else{
      #Convert to factors.
      icgc.maf$Tumor_Sample_Barcode = as.factor(x = as.character(icgc.maf$Hugo_Symbol))
      icgc.maf$Variant_Classification = as.factor(x = as.character(icgc.maf$Variant_Classification))
      icgc.maf$Variant_Type = as.factor(x = as.character(icgc.maf$Variant_Type))

      message('Summarizing..')
      maf.summary = summarizeMaf(maf = icgc.maf)

      maf.oncomat = createOncoMatrix(maf = icgc.maf)

      silent = c("3'UTR", "5'UTR", "3'Flank", "Targeted_Region", "Silent", "Intron",
                 "RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA")
      maf.silent = icgc.maf[Variant_Classification %in% silent]

      icgc.maf = MAF(data = icgc.maf, variants.per.sample = maf.summary$variants.per.sample, variant.type.summary = maf.summary$variant.type.summary,
              variant.classification.summary = maf.summary$variant.classification.summary, gene.summary = maf.summary$gene.summary,
              oncoMatrix = maf.oncomat$oncomat, numericMatrix = maf.oncomat$nummat, summary = maf.summary$summary,
              classCode = maf.oncomat$vc, maf.silent = maf.silent)
    }
  }

  return(icgc.maf)
}
