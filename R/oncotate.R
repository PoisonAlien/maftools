#' Annotates given variants using oncotator api.
#'
#' @description Takes input variants and annotates them using Broad's oncotator api (http://www.broadinstitute.org/oncotator/). Output is a dataframe of annotated variants in maf format.
#'
#' Input should be tsv file or a data.frame with first five columns in the order chr, start, end, ref_allele, alt_allele (and so on, but only first five will used, rest will be attached to resulting maf file). Note: Time consuming if input is huge.
#' Try to include necessary columns such as Tumor_Sample_Barcode along with above 5 fields. Only supports hg19/GRCh37 build.
#' @param maflite input tsv file or a data.frame with chr, start, end, ref_allele, alt_allele columns. (rest of the columns, if present will be attached to the output maf)
#' @param header logical. Whether input has a header line. Default is FALSE. Only applicable when the input is a tsv file.
#' @param basename NULL. if basename is given, annotations will be written to <basename>.maf file.
#' @return returns a data.table in maf format.
#' @examples
#' sample.var = data.frame(chromsome = c('chr4', 'chr15'), Start = c(55589774, 41961117),
#' end = c(55589774, 41961117), ref = c('A', 'TGGCTAA'), alt = c('G', '-'),
#' Tumor_Sample_Barcode = c('fake_1', 'fake2'))
#' write.table(sample.var, 'sampleVars.txt', sep='\t',quote = FALSE, row.names = FALSE)
#' ##var.maf <- oncotate(maflite = 'sampleVars.txt', header = TRUE)
#' @importFrom rjson fromJSON
#' @export

oncotate = function(maflite, header = FALSE,basename = NULL){

  #read the file
  if(is.data.frame(x = maflite)){
    m = data.table::setDT(maflite)
  }else{
    if(file.exists(maflite)){
      m = data.table::fread(maflite, stringsAsFactors = FALSE, header = header, sep='\t')
    }else{
      stop(paste0("Input file ", maflite, " doesn not exists!"))
    }
  }

  if(length(colnames(m)[colnames(m) %in% 'Tumor_Sample_Barcode']) > 0){
    Tumor_Sample_Barcode = TRUE
  }

  #paste first five columns
  anno = apply(m[,c(1:5)], 1, function(x) paste0(gsub(pattern = " ", replacement = "", x = x), collapse = "_"))

  message("Fetching annotations from Oncotator. This might take a while..")

  per_complete = as.integer(seq(1, length(anno), length.out = 4))

  pb <- txtProgressBar(min = 0, max = length(anno), style = 3)

  anno.dt = lapply(seq_along(anno), function(i){
                  setTxtProgressBar(pb, i)
                  x = anno[i]
                  rec.url = paste('http://portals.broadinstitute.org/oncotator/mutation', x, sep = '/')
                  annot = try(rjson::fromJSON(file = rec.url, unexpected.escape = "error"), silent = TRUE)
                  if(class(annot) == "try-error"){
                    annot = NA
                  }else{
                    data.table::setDT(annot, keep.rownames = TRUE)
                  }
                  annot
                })

  # failed_rows = which(sapply(anno.dt, function(x) is.na(x)))
  # if(length(failed_rows) > 0){
  #   warning(paste0("Failed to fetch annotations for ", length(failed_rows), " variants."))
  # }

  #anno.dt = anno.dt[!failed_rows]
  anno.dt = data.table::rbindlist(l = anno.dt, fill = TRUE)


  #Reformat the data according to MAF specification.
  colnames(anno.dt) = gsub(pattern = "^X",replacement = "",x = colnames(anno.dt))
  colnames(m)[1:5] = c('Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2')
  anno.dt = cbind(m, anno.dt)
  anno.dt$Center = NA
  anno.dt$Tumor_Seq_Allele1 = anno.dt$Reference_Allele
  colnames(anno.dt)[which(colnames(anno.dt) == "gene")] = "Hugo_Symbol"
  colnames(anno.dt)[which(colnames(anno.dt) == "variant_classification")] = "Variant_Classification"
  colnames(anno.dt)[which(colnames(anno.dt) == "variant_type")] = "Variant_Type"
  colnames(anno.dt)[which(colnames(anno.dt) == "HGNC_Entrez.Gene.ID.supplied.by.NCBI.")] = "Entrez_Gene_Id"
  colnames(anno.dt)[which(colnames(anno.dt) == "strand")] = "Strand"
  colnames(anno.dt)[which(colnames(anno.dt) == "build")] = "NCBI_Build"
  colnames(anno.dt)[which(colnames(anno.dt) == "strand")] = "Strand"


  col.order = c('Hugo_Symbol','Entrez_Gene_Id','Center','NCBI_Build','Chromosome','Start_Position',
                'End_Position','Strand','Variant_Classification','Variant_Type','Reference_Allele',
                'Tumor_Seq_Allele1','Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status')

  col.order = c(col.order, colnames(anno.dt)[!colnames(anno.dt) %in% col.order])

  if(length(col.order[!col.order %in% colnames(anno.dt)]) > 0){
    anno.dt = anno.dt[,col.order[col.order %in% colnames(anno.dt)], with = FALSE]
    main.cols.missing = col.order[!col.order %in% colnames(anno.dt)]
    for(i in 1:length(main.cols.missing)){
      anno.dt[,main.cols.missing[1] := NA]
    }
    anno.dt = anno.dt[,col.order, with = FALSE]
  }else{
    anno.dt = anno.dt[,col.order, with = FALSE]
  }


  if(!is.null(basename)){
    write.table(anno.df,paste(basename,'maf',sep = '.'),quote = FALSE,row.names = FALSE,sep= '\t')
  }

  return(anno.dt)
}
