#' Class MAF
#' @description S4 class for storing summarized MAF.
#' @slot data data.table of original MAF file.
#' @slot variants.per.sample table containing variants per sample
#' @slot variant.type.summary table containing variant types per sample
#' @slot variant.classification.summary table containing variant classification per sample
#' @slot gene.summary table containing variant classification per gene
#' @slot oncoMatrix character matrix of dimension n*m where n is number of genes and m is number of variants
#' @slot numericMatrix numeric matrix of dimension n*m where n is number of genes and m is number of variants
#' @slot summary table with basic MAF summary stats
#' @slot classCode mapping between numeric values in numericMatrix and Variant Classification
#' @slot maf.silent subset of main MAF containing only silent variants
#' @exportClass MAF
#' @import methods
#' @seealso \code{\link{getGeneSummary}} \code{\link{getSampleSummary}} \code{\link{getFields}}

## MAF object
MAF <- setClass(Class = 'MAF', slots =  c(data = 'data.table', variants.per.sample = 'data.table', variant.type.summary = 'data.table',
                                         variant.classification.summary = 'data.table', gene.summary = 'data.table', oncoMatrix = 'matrix',
                                         numericMatrix = 'matrix', summary = 'data.table', classCode = 'character',
                                         maf.silent = 'data.table'))

setMethod(f = 'show', signature = 'MAF', definition = function(object){
  cat(paste('An object of class ', class(object), "\n"))
  print(object@summary)
})
