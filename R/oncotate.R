#' Annotates given variants using oncotator api.
#'
#' @description Takes variants as input and annotates them using Broad's oncotator api (http://www.broadinstitute.org/oncotator/). Output is a dataframe of annotated variants in maf format.
#'
#' Input should be a five column file with chr, start, end, ref_allele, alt_allele (and so on, but only first five will used, rest will be attached to resulting maf file). Note: Time consuming if input is huge.
#' Try to include necessary columns such as Tumor_Sample_Barcode along with above 5 fields.
#' @param maflite input tsv file with chr, start, end, ref_allele, alt_allele columns. (rest of the columns, if present will be attached to the output maf)
#' @param header logical. Whether input has a header line. Default is FALSE.
#' @param basename NULL. if basename is given, annotations will be written to <basename>.maf file.
#' @return returns a dataframe in maf format.
#' @examples
#' sample.var = data.frame(chromsome = c('chr4', 'chr15'), Start = c(55589774, 41961117),
#' end = c(55589774, 41961117), ref = c('A', 'TGGCTAA'), alt = c('G', '-'),
#' Tumor_Sample_Barcode = c('fake_1', 'fake2'))
#' write.table(sample.var, 'sampleVars.txt', sep='\t',quote = FALSE, row.names = FALSE)
#' ##var.maf <- oncotate(maflite = 'sampleVars.txt', header = TRUE)
#' @importFrom rjson fromJSON
#' @export

oncotate = function(maflite, header = FALSE,basename = NULL){

  #require(package = 'rjson', quietly = TRUE)

  #create an empty data frame
  anno.df = c()

  counter = 0

  #read the file
  m = read.delim(maflite,stringsAsFactors = FALSE,header = header, sep='\t')
  #m = maflite

  if(length(colnames(m)[colnames(m) %in% 'Tumor_Sample_Barcode']) > 0){
    Tumor_Sample_Barcode = TRUE
  }

  #paste first five columns
  anno = paste(m[,1],m[,2],m[,3],m[,4],m[,5],sep = "_")

  message("fetching annotation from Oncotator. This might take a while..")
  pb = txtProgressBar(min = 0, max = length(anno), style = 3)

  for(i in 1:length(anno)){

    setTxtProgressBar(pb, value = i)
    rec = anno[i]

    #make an url for oncotator api.
    #rec.url = paste('http://www.broadinstitute.org/oncotator/mutation',rec,sep = '/') #For some reason not working
    rec.url = paste('http://portals.broadinstitute.org/oncotator/mutation',rec,sep = '/')

    #use rjason to query oncotator.
    annot = rjson::fromJSON(file = rec.url)
    anno.df = rbind(anno.df,as.data.frame(annot))
  }


  #Reformat the data according to MAF specification.
  colnames(anno.df) = gsub(pattern = "^X",replacement = "",x = colnames(anno.df))
  colnames(m)[1:5] = c('Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2')
  anno.df = cbind(m,anno.df)
  anno.df$Center = NA
  anno.df$Tumor_Seq_Allele1 = anno.df$Reference_Allele
  colnames(anno.df)[which(colnames(anno.df) == "gene")] = "Hugo_Symbol"
  colnames(anno.df)[which(colnames(anno.df) == "variant_classification")] = "Variant_Classification"
  colnames(anno.df)[which(colnames(anno.df) == "variant_type")] = "Variant_Type"
  colnames(anno.df)[which(colnames(anno.df) == "HGNC_Entrez.Gene.ID.supplied.by.NCBI.")] = "Entrez_Gene_Id"
  colnames(anno.df)[which(colnames(anno.df) == "strand")] = "Strand"
  colnames(anno.df)[which(colnames(anno.df) == "build")] = "NCBI_Build"
  colnames(anno.df)[which(colnames(anno.df) == "strand")] = "Strand"

  #get main columns
  anno.df1 = anno.df[,c('Hugo_Symbol','Entrez_Gene_Id','Center','NCBI_Build','Chromosome','Start_Position',
                        'End_Position','Strand','Variant_Classification','Variant_Type','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status')]
  #get rest of the columns
  anno.df2 = anno.df[,colnames(anno.df)[!colnames(anno.df) %in% colnames(anno.df1)]]
  #join'em
  anno.df = cbind(anno.df1,anno.df2)

  if(!is.null(basename)){
    write.table(anno.df,paste(basename,'maf',sep = '.'),quote = FALSE,row.names = FALSE,sep= '\t')
  }

  return(data.table::setDT(anno.df))
}
