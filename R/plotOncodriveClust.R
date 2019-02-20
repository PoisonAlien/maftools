#' Plots results from \code{oncodrive}
#'
#' @description Takes results from \code{oncodrive} and plots them as a scatter plot. Size of the gene shows number of clusters (hotspots), x-axis can either be an absolute number of variants
#' accumulated in these clusters or a fraction of total variants found in these clusters. y-axis is fdr values transformed into -log10 for better representation. Labels indicate Gene name with number clusters
#' observed.
#' @param  res results from \code{\link{oncodrive}}
#' @param fdrCutOff fdr cutoff to call a gene as a driver.
#' @param useFraction if TRUE uses a fraction of total variants as X-axis scale instead of absolute counts.
#' @param colCode Colors to use for indicating significant and non-signififcant genes. Default NULL
#' @param labelSize font size for labelling genes. Default 1.
#' @param bubbleSize Size for bubbles. Default 2.
#' @return a ggplot object which can be further modified.
#' @seealso \code{\link{oncodrive}}
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' laml.sig <- oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5)
#' plotOncodrive(res = laml.sig, fdrCutOff = 0.1)
#'
#' @export


plotOncodrive = function(res = NULL, fdrCutOff = 0.05, useFraction = FALSE, colCode = NULL, bubbleSize = 1, labelSize = 1){

  if(is.null(res)){
    stop('Please provide results from oncodrive.')
  }

  res$label = paste(res$Hugo_Symbol, '[',res$clusters,']', sep='')
  res$significant = ifelse(test = res$fdr < fdrCutOff, yes = 'sig', no = 'nonsig')
  res[, log_fdr := -log10(as.numeric(fdr))]

  if(is.null(colCode)){
    colCode = c('sig' = 'red', 'nonsig' = 'royalblue')
  }else{
    names(colCode)[1:2] = c('sig', 'nonsig')
  }
  res$color = ifelse(test = res$fdr < fdrCutOff, yes = colCode[1], no = colCode[2])


  par(mar = c(4, 4, 2, 2))
  if(useFraction){
    bubble_plot(plot_dat = res, x_var = "fract_muts_in_clusters", y_var = "log_fdr",
                bubble_var = "clusters", lab_dat = res[significant %in% "sig"], text_var = "label", bubble_size = bubbleSize, text_size = labelSize, col_var = "color")
    mtext(text =  "Fraction of variants within clusters", side = 1, line = 2)
  }else{
    bubble_plot(plot_dat = res, x_var = "muts_in_clusters", y_var = "log_fdr",
                bubble_var = "clusters", lab_dat = res[significant %in% "sig"], text_var = "label", bubble_size = bubbleSize, text_size = labelSize, col_var = "color")
    mtext(text =  "# of variants within clusters", side = 1, line = 2)
  }
  mtext(text =  "-log10 (fdr)", side = 2, line = 2)

}
