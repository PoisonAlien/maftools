#' Performs sample stratification based on signature contribution and performs enrichment.
#'
#' @details Performs k-means clustering to assign signature to samples and performs enrichment analysis.
#'
#' @param maf an \code{\link{MAF}} object used for signature analysis.
#' @param sig_res Signature results from \code{\link{extractSignatures}}
#' @param minMut Consider only genes with minimum this number of samples mutated. Default 5.
#' @param useCNV whether to include copy number events. Only applicable when MAF is read along with copy number data. Default TRUE if available.
#' @param fn basename for output file. Default NULL.
#' @return result list containing p-values
#' @export

signatureEnrichment = function(maf, sig_res, minMut = 5, useCNV = FALSE, fn = NULL){

  contrib = sig_res$contributions

  set.seed(seed = 1024)
  message("Running k-means for signature assignment..")
  contrib.km = kmeans(x = t(contrib), centers = nrow(contrib))
  cluster_df = as.data.frame(apply(contrib.km$centers, 2, function(x) which(x == max(x))))
  colnames(cluster_df)[1] = 'Cluster'
  #cluster_df$Cluster = paste0("Cluster_", cluster_df$Cluster)
  data.table::setDT(x = cluster_df, keep.rownames = T)
  colnames(cluster_df)[1] = 'Signature'

  xc = as.data.frame(contrib.km$cluster)
  colnames(xc)[1] = 'Cluster'
  data.table::setDT(xc, keep.rownames = TRUE)
  colnames(xc)[1] = 'Tumor_Sample_Barcode'

  xc = merge(xc, cluster_df, by = 'Cluster')
  xc[, Cluster := NULL]

  message("Performing pairwise and groupwise comparisions..")
  sig.enrich = clinicalEnrichment(maf = maf, clinicalFeature = "Signature",
                                  annotationDat = xc, minMut = minMut, useCNV = useCNV)

  sig.enrich$Signature_Assignment = xc

  cf.tbl = table(xc$Signature)
  message("Estimating mutation load per Signature..")
  mut.load = lapply(names(cf.tbl), function(x){
    tsbs = xc[Signature %in% x][,Tumor_Sample_Barcode]
    subsetMaf(maf = maf, tsb = tsbs, fields = 'Hugo_Symbol')[,.N,Tumor_Sample_Barcode]
  })

  names(mut.load) = names(cf.tbl)
  mut.load = data.table::rbindlist(mut.load, idcol = "Signature")
  mut.load[,Signature := gsub(pattern = 'Signature_', replacement = "", x = mut.load$Signature)]

  bsMax = max(unlist(lapply(split(mut.load, as.factor(mut.load$Signature)), function(x) max(boxplot.stats(x[,N])$stat))))

  mut.load.p = ggplot(data = mut.load, aes(x = Signature, y = N))+geom_boxplot(outlier.alpha = 0.6, outlier.size = 1, outlier.colour = 'gray70', fill = 'gray70')+
                      cowplot::theme_cowplot(font_size = 16, line_size = 1)+
                      theme(axis.title.x = element_text(face = "bold"), axis.text.x = element_text(face = "bold"), axis.text.y = element_text(face = "bold"), axis.title.y = element_text(face = "bold"))+
                      ylab("# mutations")+ylim(0, bsMax)

  print(mut.load.p)

  mut.load.summary = mut.load[,.(n_samples = .N, Median_mutations = median(N), Mean_mutations = mean(N)), Signature]

  sig.enrich$mutation_load = mut.load.summary

  if(!is.null(fn)){
    write.table(x = sig.enrich$pairwise_comparision, file = paste0(fn, "_pairwise_comparision.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)
    write.table(x = sig.enrich$groupwise_comparision, file = paste0(fn, "_groupwise_comparision.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)
    write.table(x = sig.enrich$Signature_Assignment, file = paste0(fn, "_Signature_Assignment.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)
    write.table(x = sig.enrich$mutation_load, file = paste0(fn, "_mutation_load.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)
    cowplot::save_plot(filename = paste0(fn, "_mutation_load.pdf"), plot = mut.load.p, base_height = 5, base_width = 4)
  }

  sig.enrich
}
