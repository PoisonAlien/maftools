#' Plot density plots from clutering results.
#' @description Plots results from inferHeterogeneity.
#' @param clusters clustering results from \code{\link{inferHeterogeneity}}
#' @param tsb sample to plot from clustering results. Default plots all samples from results.
#' @param genes genes to highlight on the plot. Can be a vector of gene names, \code{CN_altered} to label copy number altered varinats.
#'   or \code{all} to label all genes. Default NULL.
#' @param showCNvars show copy numbered altered variants on the plot. Default FALSE.
#' @param labelSize Font size for gene symbols. Default 0.8
#' @param titleSize font size for titles. Default 1.2
#' @param pointSize font size for titles. Default 1.2
#' @param savePlot If TRUE saves plot to output pdf
#' @param width plot width. Default 6.
#' @param height plot height. Default 5.
#' @param colors manual colors for clusters. Default NULL.
#' @return returns nothing.
#' @export
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' seg = system.file('extdata', 'TCGA.AB.3009.hg19.seg.txt', package = 'maftools')
#' TCGA.AB.3009.clust <- inferHeterogeneity(maf = laml, tsb = 'TCGA-AB-3009',
#' segFile = seg, vafCol = 'i_TumorVAF_WU')
#' plotClusters(TCGA.AB.3009.clust, genes = c('NF1', 'SUZ12'), showCNvars = TRUE)
#' @seealso \code{\link{inferHeterogeneity}}

plotClusters = function(clusters, tsb = NULL, genes = NULL, showCNvars = FALSE, savePlot = FALSE,
                        width = 6, height = 5, colors = NULL, titleSize = 1.2, labelSize = 0.8, pointSize = 1.2){

  clusterData = clusters$clusterData

  if(is.null(tsb)){
    tsb = as.character(unique(clusterData[,Tumor_Sample_Barcode]))
  }

  if(is.null(colors)){
    colors = RColorBrewer::brewer.pal(n = 9, name = 'Set1')
    colors = grDevices::adjustcolor(col = colors, alpha.f = 0.6)
    colors = c(colors, 'darkgray')
    names(colors) = as.character(c(as.character(1:9), 'outlier'))
  }else{
    colors = colors
  }

  for(i in 1:length(tsb)){

    tsb.dat = clusterData[Tumor_Sample_Barcode %in% tsb[i]]

    if(nrow(tsb.dat) == 0){
      stop(paste('Sample',tsb[i], 'not found'))
    }


    #CN altered and outliersregions
    tsb.dat.cn.vars = tsb.dat[cluster %in% c('CN_altered')]
    tsb.dat = tsb.dat[!cluster %in% 'CN_altered']

    tsb.dens = density(tsb.dat[!cluster %in% "outlier",t_vaf])
    nclusts = nrow(tsb.dat[,.N,cluster][!cluster %in% "outlier"])

    #lo = matrix(c(rep(1, nclusts), 2), nrow = length(c(c(rep(1, nclusts), 2))), ncol = 1, byrow = TRUE)
    #layout(mat = matrix(data = lo), heights = c(rep(1, nclusts)), 5)

    lo = layout(mat = matrix(data = c(1, 2), nrow = 2), heights = c(2, 6))

    par(mar = c(0, 3, 0, 2))
    boxplot(t_vaf ~ cluster, data = tsb.dat[!cluster %in% "outlier"], axes = FALSE,
            staplewex=0, border = colors, horizontal = TRUE, ylim = c(0, 1.2),
            lwd = 1.4, outwex=0, boxwex = 0.35)

    par(mar = c(3.5, 3, 2, 2))
    tsb.clust.spl = split(tsb.dat, as.factor(as.character(tsb.dat$cluster)))
    tsb.clust.dens = lapply(seq_along(1:nclusts), function(cd){
      cd = tsb.clust.spl[[cd]]
      cd.dens = density(cd[,t_vaf])
      cd.dens
    })
    yl = c(0, round(max(unlist(lapply(tsb.clust.dens, function(x) max(x$y, na.rm = TRUE)))), digits = 2))

    plot(0, 0, xlim = c(0, 1.2), ylim = yl,
         main = "", lwd = 2, xlab = "", ylab = "", axes = FALSE, pch = NA)
    points(tsb.dens, type = "l", lwd = 1.5,
           col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))


    for(clust in 1:nclusts){
      points(tsb.clust.dens[[clust]], type = "l", lwd = 3, col = colors[clust])
    }

    axis(side = 2, at = yl, las = 2, hadj = 0.85, lwd = 3, font.axis = 2, cex = 1.5, font = 2)
    axis(side = 1, at = seq(0, 1, by = 0.25), padj = -0.5, lwd = 3, font.axis = 2, cex = 1.5, font = 2)
    points(x = tsb.dat[,t_vaf], y = rep(0, nrow(tsb.dat)), col = colors[tsb.dat$cluster], pch = 16, cex = pointSize)
    title(main = tsb[i], font = 2, adj = 0, line = 1, cex.main = titleSize)
    title(main = paste0("MATH: ", round(unique(tsb.dat[,MATH]), digits = 2)),
          font.main = 4, adj = 1, line = 1, cex.main = 0.8*titleSize)

    legend(x = "topright", legend = names(colors[unique(tsb.dat[,cluster])]), col = colors[unique(tsb.dat[,cluster])],
           border = NA, bty = "n", pch = 16, xpd = TRUE, ncol = 1,
           cex = 1.2, pt.cex = 1.5)

    mtext(text = "VAF", side = 1, line = 2, font = 2, cex = 1)

    if(showCNvars){
      if(nrow(tsb.dat.cn.vars) > 0){
        points(x = tsb.dat.cn.vars[,t_vaf], y = rep(0, nrow(tsb.dat.cn.vars)), col = "#00000099", pch = 16)
      }
    }

    #Are there genes to highlight?
    if(!is.null(genes)){
      if(genes[1] == 'CN_altered'){
        if(nrow(tsb.dat.cn.vars) > 0){
          ypos = seq(0.15*max(yl, na.rm = TRUE), 0.3*max(yl, na.rm = TRUE), length.out = nrow(tsb.dat.cn.vars))
          segments(x0 = tsb.dat.cn.vars[,t_vaf], y0 = 0, x1 = tsb.dat.cn.vars[,t_vaf], y1 = ypos, lwd = 2, col = "gray70")
          text(x = tsb.dat.cn.vars[,t_vaf], y = ypos, labels = tsb.dat.cn.vars[,Hugo_Symbol], font = 2, cex = labelSize)
        }
      }else if(length(genes) > 0){
        ypos = seq(0.15*max(yl, na.rm = TRUE), 0.3*max(yl, na.rm = TRUE), length.out = length(genes))
        genesDat = rbind(tsb.dat, tsb.dat.cn.vars, fill = TRUE)[Hugo_Symbol %in% genes]
        if(nrow(genesDat) > 0){
          segments(x0 = genesDat[,t_vaf], y0 = 0, x1 = genesDat[,t_vaf], y1 = ypos, lwd = 2, col = "gray70")
          text(x = genesDat[,t_vaf], y = ypos, labels = genesDat[,Hugo_Symbol], font = 2, cex = labelSize)
        }
      }
    }
  }
}
