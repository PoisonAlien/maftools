#' Class MAF
#' @description S4 class for storing summarized MAF.
#' @slot data data.table of MAF file containing all non-synonymous variants.
#' @slot variants.per.sample table containing variants per sample
#' @slot variant.type.summary table containing variant types per sample
#' @slot variant.classification.summary table containing variant classification per sample
#' @slot gene.summary table containing variant classification per gene
#' @slot summary table with basic MAF summary stats
#' @slot maf.silent subset of main MAF containing only silent variants
#' @slot clinical.data clinical data associated with each sample/Tumor_Sample_Barcode in MAF.
#' @exportClass MAF
#' @import methods
#' @seealso \code{\link{getGeneSummary}} \code{\link{getSampleSummary}} \code{\link{getFields}}

## MAF object
setClass(Class = 'MAF', slots =  c(data = 'data.table', variants.per.sample = 'data.table', variant.type.summary = 'data.table',
                                         variant.classification.summary = 'data.table', gene.summary = 'data.table',
                                         summary = 'data.table', maf.silent = 'data.table', clinical.data = 'data.table'))

setMethod(f = 'show', signature = 'MAF', definition = function(object){
  cat(paste('An object of class ', class(object), "\n"))
  print(object@summary)
})


#' Construct an MAF object
#' @description Constructor function which takes non-synonymous, and synonymous variants along with an optional clinical information and generates an MAF object
#' @param nonSyn non-synonymous variants as a data.table or any object that can be coerced into a data.table (e.g: data.frame, GRanges)
#' @param syn synonymous variants as a data.table or any object that can be coerced into a data.table (e.g: data.frame, GRanges)
#' @param clinicalData Clinical data associated with each sample/Tumor_Sample_Barcode in MAF. Could be a text file or a data.frame. Requires at least a column with the name `Tumor_Sample_Barcode` Default NULL.
#' @param verbose Default TRUE
#' @export
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml_dt = data.table::fread(input = laml.maf)
#' laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') #Clinical data
#'# Just for demonstration
#' nsyn_vars = laml_dt[Variant_Classification %in% "Missense_Mutation"]
#' syn_vars = laml_dt[Variant_Classification %in% "Silent"]
#' maftools::MAF(nonSyn = nsyn_vars, syn = syn_vars, clinicalData = laml.clin)

MAF = function(nonSyn = NULL, syn = NULL, clinicalData = NULL, verbose = TRUE) {

  if(is.null(nonSyn)){
    nonSyn = data.table::data.table()
  }else if(!is(object = nonSyn, class2 = "data.frame")){
    nonSyn = data.table::as.data.table(nonSyn)
  }else{
    data.table::setDT(nonSyn)
  }

  if(is.null(syn)){
    syn = data.table::data.table()
  }else if(!is(object = syn, class2 = "data.frame")){
    syn = data.table::as.data.table(syn)
  }else{
    data.table::setDT(syn)
  }

  mafSummary = summarizeMaf(maf = nonSyn, anno = clinicalData, chatty = verbose)

  new(Class = "MAF", data = nonSyn, maf.silent = syn, clinical.data = mafSummary$sample.anno, variants.per.sample = mafSummary$variants.per.sample, variant.type.summary = mafSummary$variant.type.summary,
      variant.classification.summary = mafSummary$variant.classification.summary, gene.summary = mafSummary$gene.summary,
      summary = mafSummary$summary)

}


#
#' Convert MAF to MultiAssayExperiment object
#' @description Generates an object of class \code{MultiAssayExperiment} from \code{MAF} object
#' @param m an \code{MAF} object
#' @export
#' @examples
#' laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
#' laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')
#' laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
#' maf2mae(laml)

maf2mae = function(m = NULL){

  if(any(!requireNamespace("RaggedExperiment") | !requireNamespace("MultiAssayExperiment"))){
    stop("Converting to MultiAssayExperiment requires RaggedExperiment and MultiAssayExperiment packages!")
  }

  if(is.null(m)){
    stop("Missing input MAF")
  }

  # Convert clinical data to coldata
  maf_coldata = S4Vectors::DataFrame(getClinicalData(m))
  rownames(maf_coldata) = maf_coldata$Tumor_Sample_Barcode

  # Make RaggedExperiment for non-synonymous variants
  maf_nsyn_gr = GenomicRanges::GRanges(
    seqnames = m@data$Chromosome,
    ranges = IRanges::IRanges(
      start = m@data$Start_Position,
      end = m@data$End_Position
    ), strand = "*",
    S4Vectors::DataFrame(m@data[, setdiff(x = colnames(m@data),
                                          y = c("Chromosome", "Start_Position", "End_Position")), with = FALSE])
  )
  maf_nsyn_gr = split(maf_nsyn_gr, as.factor(maf_nsyn_gr$Tumor_Sample_Barcode))
  maf_nsyn_re = RaggedExperiment::RaggedExperiment(maf_nsyn_gr, colData = maf_coldata)

  # Make RaggedExperiment for synonymous variants
  maf_syn_gr = GenomicRanges::GRanges(
    seqnames = m@maf.silent$Chromosome,
    ranges = IRanges::IRanges(
      start = m@maf.silent$Start_Position,
      end = m@maf.silent$End_Position
    ), strand = "*",
    S4Vectors::DataFrame(m@maf.silent[, setdiff(x = colnames(m@maf.silent),
                                                y = c("Chromosome", "Start_Position", "End_Position")), with = FALSE])
  )
  maf_syn_gr = split(maf_syn_gr, as.factor(maf_syn_gr$Tumor_Sample_Barcode))
  maf_syn_re = RaggedExperiment::RaggedExperiment(maf_syn_gr, colData = maf_coldata)


  # Create MultiAssayExperiment object with metadata
  MultiAssayExperiment::MultiAssayExperiment(
    experiments = list("maf_nonSyn" = maf_nsyn_re, "maf_syn" = maf_syn_re),
    metadata = list(
      summary = m@summary,
      variant.classification.summary = m@variant.classification.summary,
      variant.type.summary = m@variant.type.summary,
      variants.per.sample = m@variants.per.sample,
      gene.summary = m@gene.summary
    ),
    colData = maf_coldata
  )
}
