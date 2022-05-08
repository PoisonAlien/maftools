#' Plot density plots from clutering results.
#' @description Plots results from inferHeterogeneity.
#' @param clusters clustering results from \code{\link{inferHeterogeneity}}
#' @param tsb sample to plot from clustering results. Default plots all samples from results.
#' @param genes genes to highlight on the plot. Can be a vector of gene names, \code{CN_altered} to label copy number altered varinats.
#'   or \code{all} to label all genes. Default NULL.
#' @param showCNvars show copy numbered altered variants on the plot. Default FALSE.
#' @param colors manual colors for clusters. Default NULL.
#' @return returns nothing.
#' @export
#' @examples
#' \dontrun{
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' seg = system.file('extdata', 'TCGA.AB.3009.hg19.seg.txt', package = 'maftools')
#' TCGA.AB.3009.clust <- inferHeterogeneity(maf = laml, tsb = 'TCGA-AB-3009',
#' segFile = seg, vafCol = 'i_TumorVAF_WU')
#' plotClusters(TCGA.AB.3009.clust, genes = c('NF1', 'SUZ12'), showCNvars = TRUE)
#' }
#' @seealso \code{\link{inferHeterogeneity}}

plotClusters = function(clusters, tsb = NULL, genes = NULL, showCNvars = FALSE, colors = NULL){

  clusterData = clusters$clusterData

  if(is.null(tsb)){
    tsb = as.character(unique(clusterData[,Tumor_Sample_Barcode]))
  }

  if(is.null(colors)){
    colors = RColorBrewer::brewer.pal(n = 9, name = 'Set1')
    colors = c(colors, 'darkgray', 'black')
    names(colors) = as.character(c(as.character(1:9), 'outlier', 'CN_altered'))
  }else{
    colors = colors
  }

  for(i in 1:length(tsb)){

    tsb.dat = clusterData[Tumor_Sample_Barcode %in% tsb[i]]

    if(nrow(tsb.dat) == 0){
      stop(paste('Sample',tsb[i], 'not found'))
    }


    #CN altered and outliersregions
    tsb.dat.cn.vars = tsb.dat[cluster == c('CN_altered')]
    #tsb.dat = tsb.dat[!cluster == 'CN_altered']
    tsb_nclust = as.numeric(unique(tsb.dat[!cluster %in% c("outlier", "CN_altered"),cluster]))
    lo_mat = matrix(c(tsb_nclust, max(tsb_nclust)+1) , ncol = 1)
    #return(lo_mat)
    graphics::layout(mat = lo_mat, heights = c(rep(0.5, length(tsb_nclust)), 4))

    tsb.dat.spl = split(tsb.dat, as.factor(as.character(tsb.dat$cluster)))
    for(cl in seq_along(tsb_nclust)){
      par(mar = c(0, 3, 0, 2))
      boxplot(tsb.dat.spl[[cl]][,t_vaf], axes = FALSE, ylim = c(0, 1), horizontal = TRUE,
              xaxt="n", boxwex=0.75, outline=FALSE, lty=1, lwd = 1.4, outwex=0, staplewex=0,
              border = colors[cl])
    }

    tsb_dens = density(x = tsb.dat[!cluster %in% c("outlier", "CN_altered"), t_vaf])
    y_ax_ticks = pretty(x = range(tsb_dens$y))

    par(mar = c(4, 3, 3, 2))
    plot(NA, NA, axes = FALSE, xlab = NA, ylab = NA, xlim = c(0, 1),
         ylim = y_ax_ticks[c(1, length(y_ax_ticks))])
    lines(tsb_dens, lwd = 1.5)
    abline(h = y_ax_ticks, v = seq(0, 1, 0.2), lty = 2, col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.5))
    axis(side = 1, at = seq(0, 1, 0.2), lwd = 1.5, cex = 1)
    axis(side = 2, at = y_ax_ticks, lwd = 1.5, las = 2, cex = 1)
    points(y = rep(0, length(tsb.dat[, t_vaf])),
           x = tsb.dat[, t_vaf], pch = 16, col = colors[tsb.dat[, cluster]], cex = 1)
    title(main = tsb[i], adj = 0, line = 2)
    title(main = paste0("MATH: ", round(unique(tsb.dat[!is.na(MATH)][,MATH]), digits = 3)), cex.main = 1,
          font.main= 3, col.main= "black", cex.sub = 1, font.sub = 3, line = 0.75, adj = 0)

    #Are there genes to highlight?
    if(!is.null(genes)){
      if("CN_altered" %in% genes){
        if(nrow(tsb.dat.cn.vars) > 0){
          segments(x0 = tsb.dat.cn.vars[, t_vaf], y0 = 0, x1 = tsb.dat.cn.vars[, t_vaf], y1 = 0.1 * max(tsb_dens$y, na.rm = TRUE), lwd = 1)
          text(x = tsb.dat.cn.vars[, t_vaf], y = 0.1 * max(tsb_dens$y, na.rm = TRUE), srt = 90, font = 3, labels = tsb.dat.cn.vars[, Hugo_Symbol], adj = 0)
        }
      }else if('all' %in% genes){
        segments(x0 = tsb.dat[, t_vaf], y0 = 0, x1 = tsb.dat[, t_vaf], y1 = 0.1 * max(tsb_dens$y, na.rm = TRUE), lwd = 1)
        text(x = tsb.dat[, t_vaf], y = 0.1 * max(tsb_dens$y, na.rm = TRUE),
             srt = 90, font = 3, labels = tsb.dat[, Hugo_Symbol], adj = 0, cex = 1)
      }else{
        genesDat = tsb.dat[Hugo_Symbol %in% genes]
        if(nrow(genesDat) > 0){
          segments(x0 = genesDat[, t_vaf], y0 = 0, x1 = genesDat[, t_vaf], y1 = 0.1 * max(tsb_dens$y, na.rm = TRUE), lwd = 1)
          text(x = genesDat[, t_vaf], y = 0.1 * max(tsb_dens$y, na.rm = TRUE), srt = 90, font = 3, labels = genesDat[, Hugo_Symbol], adj = 0)
      }
      }
    }

    legend(x = "topright", legend = unique(names(colors[tsb.dat[, cluster]])), bty = "n",
           col = unique(colors[tsb.dat[, cluster]]), pch = 16, border = NA, cex = 1,
           title = "Cluster")
    mtext(text = "VAF", side = 1, cex = 1.2, line = 2.2)
  }
}
