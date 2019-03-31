#' Performs sample stratification based on signature contribution and enrichment analysis.
#'
#' @description Performs k-means clustering to assign signature to samples and performs enrichment analysis.
#'
#' @param maf an \code{\link{MAF}} object used for signature analysis.
#' @param sig_res Signature results from \code{\link{extractSignatures}}
#' @param minMut Consider only genes with minimum this number of samples mutated. Default 5.
#' @param useCNV whether to include copy number events. Only applicable when MAF is read along with copy number data. Default TRUE if available.
#' @param fn basename for output file. Default NULL.
#' @return result list containing p-values
#' @export
#' @seealso \code{\link{plotEnrichmentResults}}

signatureEnrichment = function(maf, sig_res, minMut = 5, useCNV = FALSE, fn = NULL){

  contrib = sig_res$contributions
  contrib = t(contrib)

  set.seed(seed = 1024)
  message("Running k-means for signature assignment..")
  contrib.km = kmeans(x = contrib, centers = ncol(contrib))
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

  sig.mean.stat = lapply(unique(as.character(xc[,Signature])), function(x){
              data.frame(sig_mean =
                apply(contrib[xc[Signature %in% x, Tumor_Sample_Barcode],, drop = FALSE], 2, mean)
              )
            })

  sig.sd.stat = lapply(unique(as.character(xc[,Signature])), function(x){
    data.frame(sig_sd =
      apply(contrib[xc[Signature %in% x, Tumor_Sample_Barcode],, drop = FALSE], 2, sd)/sqrt(nrow(xc[Signature %in% x]))
    )
  })

  names(sig.mean.stat) = names(sig.sd.stat)  = unique(as.character(xc[,Signature]))

  sig.mean.stat = as.data.frame(sig.mean.stat)
  colnames(sig.mean.stat) = unique(as.character(xc[,Signature]))

  sig.sd.stat = as.data.frame(sig.sd.stat)
  colnames(sig.sd.stat) = unique(as.character(xc[,Signature]))


  message("Performing pairwise and groupwise comparisions..")
  sig.enrich = clinicalEnrichment(maf = maf, clinicalFeature = "Signature",
                                  annotationDat = xc, minMut = minMut, useCNV = useCNV)

  sig.enrich$Signature_Assignment = xc

  cf.tbl = table(xc$Signature)
  message("Estimating mutation load and signature exposures..")
  mut.load = lapply(names(cf.tbl), function(x){
    tsbs = xc[Signature %in% x][,Tumor_Sample_Barcode]
    subsetMaf(maf = maf, tsb = tsbs, fields = 'Hugo_Symbol', mafObj = FALSE)[,.N,Tumor_Sample_Barcode]
  })

  names(mut.load) = names(cf.tbl)
  mut.load = data.table::rbindlist(mut.load, idcol = "Signature")
  mut.load[,Signature := gsub(pattern = 'Signature_', replacement = "", x = mut.load$Signature)]

  bsMax = max(unlist(lapply(split(mut.load, as.factor(mut.load$Signature)), function(x) max(boxplot.stats(x[,N])$stat))))

  cols = RColorBrewer::brewer.pal(n = 9, name = "Pastel1")[1:ncol(sig.mean.stat)]
  names(cols) = rownames(sig.mean.stat)

  add_legend <- function(...) {
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
                mar=c(0, 0, 0, 0), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend(...)
  }


  par(bty="n", mgp = c(0.5,0.5,0), las=1, tcl=-.25, font.main=4, xpd=TRUE, mar = c(3.4,2,2,4), mfrow = c(1, 2))
  #layout(matrix(1:1, ncol=2, byrow = TRUE), heights=c(1,1,1,0.2))
  title_size = 1
  par(mar = c(3.5,2.5,2,2))
  b = barplot(as.matrix(sig.mean.stat), ylim = c(0, 1 + max(sig.sd.stat, na.rm = TRUE)), col = cols,
              axes = FALSE, border = 0.1, xaxt = "n")
  for(i in 1:ncol(sig.sd.stat)){
    segments(x0 = b[i], y0 = cumsum(sig.mean.stat[,i]) - sig.sd.stat[,i],
             x1 = b[i], y1 = cumsum(sig.mean.stat[,i]) + sig.sd.stat[,i],
             lwd = 1.5)
  }
  axis(side = 2, at = seq(0, 1, 0.25), labels = seq(0, 1, 0.25),
       lwd = 1.2, font.axis = 1, cex = 1.5, font = 1)
  mtext(text = gsub(pattern = "Signature_", replacement = "", x = colnames(sig.mean.stat)),
        side = 1, at = b, font  = 1)
  mtext(text = "k-mean signature cluster", side = 1, line = 1.5, font = 1, cex = 0.9)
  title(main = 'Avg. signature exposure', cex.main = title_size, adj = 0)

  par(mar = c(3.4,2.5,2,2))
  boxplot(at = 1:nrow(sig.mean.stat), N ~ Signature, data = mut.load,  xaxt="n", col = cols,
          boxwex=0.6, outline = FALSE, lty=1, outwex=0, staplewex=0, axes = FALSE, ylim = c(0, bsMax), xlab = NA, ylab = NA)
  axis(side = 2, at = as.integer(seq(0, bsMax, length.out = 5)), lwd = 1, font = 1)
  mtext(text = c(gsub(pattern = "Signature_", replacement = "", x = names(cols))),
                 side = 1, at = 1:nrow(sig.mean.stat), font  = 1)
  mtext(text = c("N:", as.numeric(cf.tbl)),
        side = 1, line = 1.5, font = 1, at = 0:nrow(sig.mean.stat))
  #mtext(text = "Signature", side = 1, line = 2, font = 2)

  title(main = 'Mutation load', cex.main = title_size, adj = 0)

  # add_legend("topright", pt.lwd = 2,
  #            legend = names(cols), fill = cols,
  #            bty = "n", cex = 1, border=NA, xpd = TRUE, text.font = 2)

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
  }

  sig.enrich
}
