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
MAF <- setClass(Class = 'MAF', slots =  c(data = 'data.table', variants.per.sample = 'data.table', variant.type.summary = 'data.table',
                                         variant.classification.summary = 'data.table', gene.summary = 'data.table',
                                         summary = 'data.table', maf.silent = 'data.table', clinical.data = 'data.table'))

setMethod(f = 'show', signature = 'MAF', definition = function(object){
  cat(paste('An object of class ', class(object), "\n"))
  print(object@summary)
})

create_maf = function(maf_summary){
  req_fields = c("data", "variants.per.sample", "variant.type.summary", "variant.classification.summary",
                 "gene.summary", "summary", "maf.silent", "clinical.data")

  if(!all(names(mafSummary) %in% req_fields)){
    missing_fields = setdiff(x = req_fields, names(mafSummary))
    stop("Missing required fields:\n", paste(missing_fields, collapse = "\n"))
  }

  maf_summary = lapply(maf_summary, function(x){
    if(!is(object = x, class2 = "data.frame")){
      #try to coerce: (will lead to error if fails)
      x = data.table::as.data.table(x)
    }
    x
  })

  do.call(maftools:::MAF, maf_summary)
}
