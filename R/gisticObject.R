#' Class GISTIC
#' @description S4 class for storing summarized MAF.
#' @slot data data.table of summarized GISTIC file.
#' @slot cnv.summary table containing alterations per sample
#' @slot cytoband.summary table containing alterations per cytoband
#' @slot gene.summary table containing alterations per gene
#' @slot cnMatrix character matrix of dimension n*m where n is number of genes and m is number of samples
#' @slot numericMatrix numeric matrix of dimension n*m where n is number of genes and m is number of samples
#' @slot summary table with basic GISTIC summary stats
#' @slot classCode mapping between numeric values in numericMatrix and copy number events.
#' @exportClass GISTIC
#' @import methods
#' @seealso \code{\link{getGeneSummary}} \code{\link{getSampleSummary}} \code{\link{getCytobandSummary}}

## MAF object
GISTIC <- setClass(Class = 'GISTIC', slots =  c(data = 'data.table', cnv.summary = 'data.table',
                                             cytoband.summary = 'data.table', gene.summary = 'data.table', cnMatrix = 'matrix',
                                          numericMatrix = 'matrix', summary = 'data.table', classCode = 'character'
                                          ))

setMethod(f = 'show', signature = 'GISTIC', definition = function(object){
  cat(paste('An object of class ', class(object), "\n"))
  print(object@summary)
})
