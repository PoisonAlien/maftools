#' Draw two oncoplots side by side for cohort comparision.
#' @details Draws two oncoplots side by side to display difference between two cohorts.
#'
#' @param m1 first \code{\link{MAF}} object
#' @param m2 second \code{\link{MAF}} object
#' @param genes draw these genes. Default plots top 5 mutated genes from two cohorts.
#' @param colors named vector of colors for each Variant_Classification.
#' @param removeNonMutated Logical. If \code{TRUE} removes samples with no mutations in the oncoplot for better visualization. Default \code{FALSE}.
#' @param m1Name optional name for first cohort
#' @param m2Name optional name for second cohort
#' @export
#' @examples
#' #' ##Primary and Relapse APL
#' primary.apl <- system.file("extdata", "APL_primary.maf.gz", package = "maftools")
#' relapse.apl <- system.file("extdata", "APL_relapse.maf.gz", package = "maftools")
#' ##Read mafs
#' primary.apl <- read.maf(maf = primary.apl)
#' relapse.apl <- read.maf(maf = relapse.apl)
#' ##Plot
#' coOncoplot(m1 = primary.apl, m2 = relapse.apl, m1Name = 'Primary APL', m2Name = 'Relapse APL')
#' dev.off()
#' @return Returns nothing. Just draws plot.

coOncoplot = function(m1, m2, genes = NULL, colors = NULL, removeNonMutated = FALSE, m1Name = NULL, m2Name = NULL){

  if(is.null(genes)){
    m1.genes = getGeneSummary(m1)[1:5]
    m2.genes = getGeneSummary(m2)[1:5]
    mdt = merge(m1.genes[,.(Hugo_Symbol, MutatedSamples)], m2.genes[,.(Hugo_Symbol, MutatedSamples)], by = 'Hugo_Symbol', all = TRUE)
    mdt$MutatedSamples.x[is.na(mdt$MutatedSamples.x)] = 0
    mdt$MutatedSamples.y[is.na(mdt$MutatedSamples.y)] = 0
    mdt$max = apply(mdt[,.(MutatedSamples.x, MutatedSamples.y)], 1, max)
    mdt = mdt[order(max, decreasing = TRUE)]

    genes = mdt[,Hugo_Symbol]
  }

  m1.sampleSize = m1@summary[3, summary]
  m2.sampleSize = m2@summary[3, summary]


  if(is.null(m1Name)){
    m1Name = 'M1'
  }

  m1Name = paste(m1Name, ' (n = ' , m1.sampleSize, ')',sep = '')

  if(is.null(m2Name)){
    m2Name = 'M2'
  }

  m2Name = paste(m2Name, ' (n = ' , m2.sampleSize, ')',sep = '')

  m1.oc = getOncoPlot(maf = m1, genes = genes, removeNonMutated = removeNonMutated,
                      colors = colors, showGenes = TRUE, left = TRUE, hmName = m1Name)
  m2.oc = getOncoPlot(maf = m2, genes = genes, removeNonMutated = removeNonMutated,
                      colors = colors, showGenes = FALSE, left = FALSE, hmName = m2Name)

  oc.list = m1.oc[[1]] + m2.oc[[1]]

  #ComplexHeatmap::draw(oc.list)

  tn = unlist(unique(c(m1.oc[[2]], m2.oc[[2]])))
  tn = unique(unlist(strsplit(x = tn, split = ';', fixed = TRUE)))
  tc = c(m1.oc[[3]], m2.oc[[3]])

  legend = grid::legendGrob(labels = tn,  pch = 15, gp = grid::gpar(col = tc[tn]), nrow = 2)

  suppressWarnings( ComplexHeatmap::draw(oc.list, newpage = FALSE, annotation_legend_side = "bottom", annotation_legend_list = list(legend)) )

}
