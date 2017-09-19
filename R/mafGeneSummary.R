#' Summary statistics of MAF
#' @description Summarizes genes and samples irrespective of the type of alteration. This is different from \code{\link{getSampleSummary}} and \code{\link{getGeneSummary}} which returns summaries of only non-synonymous variants.
#' @details This function takes MAF object as input and returns summary table.
#' @return Returns a list of summarized tables
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' mafSummary(maf = laml)
#'
#' @seealso \code{\link{getGeneSummary}} \code{\link{getSampleSummary}}
#' @export

mafSummary = function(maf){
  x = summarizeMaf(maf = subsetMaf(maf = maf, fields = 'Hugo_Symbol', includeSyn = TRUE), chatty = FALSE)[c("variant.classification.summary", "gene.summary", "variants.per.sample", "variant.type.summary")]
  return(x)
}
