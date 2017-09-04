#' Class MAF
#' @description S4 class for storing summarized MAF.
#' @slot data data.table of MAF file containing all non-synonymous variants.
#' @slot variants.per.sample table containing variants per sample
#' @slot variant.type.summary table containing variant types per sample
#' @slot variant.classification.summary table containing variant classification per sample
#' @slot gene.summary table containing variant classification per gene
#' @slot summary table with basic MAF summary stats
#' @slot maf.silent subset of main MAF containing only silent variants
#' @slot sample.anno annotations associated with each sample/Tumor_Sample_Barcode in MAF. Ideally annotation would contain clinical data, survival information and other features features associated with samples.
#' @exportClass MAF
#' @import methods
#' @seealso \code{\link{getGeneSummary}} \code{\link{getSampleSummary}} \code{\link{getFields}}

## MAF object
MAF <- setClass(Class = 'MAF', slots =  c(data = 'data.table', variants.per.sample = 'data.table', variant.type.summary = 'data.table',
                                         variant.classification.summary = 'data.table', gene.summary = 'data.table',
                                         summary = 'data.table', maf.silent = 'data.table', sample.anno = 'data.table'))

setMethod(f = 'show', signature = 'MAF', definition = function(object){
  cat(paste('An object of class ', class(object), "\n"))
  print(object@summary)
})
