#' Plots results from \code{oncodrive}
#'
#' @description Takes results from \code{oncodrive} and plots them as a scatter plot. Size of the gene shows number of clusters (hotspots), x-axis can either be an absolute number of variants
#' accumulated in these clusters or a fraction of total variants found in these clusters. y-axis is fdr values transformed into -log10 for better representation. Labels indicate Gene name with number clusters
#' observed.
#' @param  res results from \code{\link{oncodrive}}
#' @param fdrCutOff fdr cutoff to call a gene as a driver.
#' @param useFraction if TRUE uses a fraction of total variants as X-axis scale instead of absolute counts.
#' @return a ggplot object which can be further modified.
#' @seealso \code{\link{oncodrive}}
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf, removeSilent = TRUE, useAll = FALSE)
#' laml.sig <- oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5)
#' plotOncodrive(res = laml.sig, fdrCutOff = 0.1)
#'
#' @export


plotOncodrive = function(res = NULL, fdrCutOff = 0.05, useFraction = FALSE){

  if(is.null(res)){
    stop('Please provide results from oncodrive.')
  }


  res$label = paste(res$Hugo_Symbol, '[',res$clusters,']', sep='')
  res$significant = ifelse(test = res$fdr < fdrCutOff, yes = 'sig', no = 'nonsig')

  if(useFraction){
    p = ggplot(data = res, aes(x = fract_muts_in_clusters, y = -log10(fdr), size = clusters, color = significant))+
      geom_point(alpha = 0.9)+cowplot::theme_cowplot()+theme(legend.position = 'NONE')+scale_color_manual(values = c('sig' = 'maroon', 'nonsig' = 'blue'))+
      ggrepel::geom_text_repel(data = res[fdr < fdrCutOff], aes(x = fract_muts_in_clusters, y = -log10(fdr), label = label, size = 2))+
      xlab('Fraction of mutations in clusters')+cowplot::background_grid(major = 'xy')
  }else{
    p = ggplot(data = res, aes(x = muts_in_clusters, y = -log10(fdr), size = clusters, color = significant))+
      geom_point()+cowplot::theme_cowplot()+theme(legend.position = 'NONE')+scale_color_manual(values = c('sig' = 'maroon', 'nonsig' = 'blue'))+
      ggrepel::geom_text_repel(data = res[fdr < fdrCutOff], aes(x = muts_in_clusters, y = -log10(fdr), label = label, size = 2))+
      xlab('Number of mutations in clusters')+cowplot::background_grid(major = 'xy')
  }

  print(p)
  #return(p)
}
