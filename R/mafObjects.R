#' @import methods
#' @exportMethod getGeneSummary
#' @exportMethod getSampleSummary
#' @export

## MAF object
MAF <- setClass(Class = 'MAF', slots =  c(data = 'data.table', variants.per.sample = 'data.table', variant.type.summary = 'data.table',
                                         variant.classification.summary = 'data.table', gene.summary = 'data.table', oncoMatrix = 'matrix',
                                         numericMatrix = 'matrix', summary = 'data.table', classCode = 'character',
                                         maf.silent = 'data.table'))

## Accessor methods
setGeneric(name = "getGeneSummary", function(x, ...) standardGeneric("getGeneSummary"))
setMethod(f = "getGeneSummary",signature = "MAF", function(x) x@gene.summary)

setGeneric(name = "getSampleSummary", function(x, ...) standardGeneric("getSampleSummary"))
setMethod(f = "getSampleSummary",signature = "MAF", function(x) x@variant.classification.summary)

setMethod(f = 'show', signature = 'MAF', definition = function(object){
  cat(paste('An object of class ', class(object), "\n"))
  print(object@summary)
})
