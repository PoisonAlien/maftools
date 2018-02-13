#' Plot density plots from clutering results.
#' @description Plots results from inferHeterogeneity.
#' @param clusters clustering results from \code{\link{inferHeterogeneity}}
#' @param tsb sample to plot from clustering results. Default plots all samples from results.
#' @param genes genes to highlight on the plot. Can be a vector of gene names, \code{CN_altered} to label copy number altered varinats.
#'   or \code{all} to label all genes. Default NULL.
#' @param showCNvars show copy numbered altered variants on the plot. Default FALSE.
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

plotClusters = function(clusters, tsb = NULL, genes = NULL, showCNvars = FALSE, savePlot = FALSE, width = 6, height = 5, colors = NULL){

  clusterData = clusters$clusterData

  if(is.null(tsb)){
    tsb = as.character(unique(clusterData[,Tumor_Sample_Barcode]))
  }

  if(is.null(colors)){
    colors = RColorBrewer::brewer.pal(n = 9, name = 'Set1')
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
    tsb.dat.cn.vars = tsb.dat[cluster == c('CN_altered')]
    tsb.dat = tsb.dat[!cluster == 'CN_altered']
    #Top density plot
    tsb.dat.dens = ggplot(tsb.dat, aes(t_vaf))+geom_density(data = tsb.dat[cluster != 'outlier'], size = 1, alpha = 0.3)+geom_point(aes(y = 0, x = t_vaf, color = cluster), size = 3, alpha = 0.6)+
      cowplot::theme_cowplot(line_size = 1)+theme(legend.position = 'bottom', plot.title = element_text(size = 12), axis.text.x = element_text(face = "bold"), axis.text.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"), axis.title.y = element_blank())+xlab('VAF')+xlim(0,1)+ylab('')+
      cowplot::background_grid('xy')+
      scale_colour_manual(values = colors)+
      ggtitle(tsb[i], subtitle = paste0("MATH: ", round(unique(tsb.dat[,MATH]), digits = 3)))

    #If any copy number altered variants, mark them with dark dots.
    if(nrow(tsb.dat.cn.vars) > 0){
      if(showCNvars){
        tsb.dat.dens = tsb.dat.dens+geom_point(data = tsb.dat.cn.vars, aes(y = 0, x = t_vaf), size = 3,alpha = 0.6)
      }
    }

    #Are there genes to highlight?
    if(!is.null(genes)){
      if(length(genes) > 1){
        genesDat = rbind(tsb.dat, tsb.dat.cn.vars, fill = TRUE)
        genesDat = genesDat[Hugo_Symbol %in% genes]
        #genesDat = dplyr::filter(.data = rbind(tsb.dat, tsb.dat.cn.vars, fill = TRUE), filter = Hugo_Symbol %in% genes)
        if(nrow(genesDat) > 0){
          tsb.dat.dens = tsb.dat.dens+ggrepel::geom_text_repel(data = genesDat,
                                                               aes(label = Hugo_Symbol, x = t_vaf, y = 0), force = 10, nudge_y = 0.5)
        }
      }else if(genes == 'all'){
        tsb.dat.dens = tsb.dat.dens+ggrepel::geom_text_repel(data = rbind(tsb.dat, tsb.dat.cn.vars, fill = TRUE),
                                                             aes(label = Hugo_Symbol, x = t_vaf, y = 0), force = 10, nudge_y = 0.5)
      }else if(genes == 'CN_altered'){
        if(nrow(tsb.dat.cn.vars) > 0){
          tsb.dat.dens = tsb.dat.dens+ggrepel::geom_text_repel(data = tsb.dat.cn.vars,
                                                               aes(label = Hugo_Symbol, x = t_vaf, y = 0), force = 10, nudge_y = 0.5)
        }
      }
    }

    #Top boxplot
    tsb.dat.bar = ggplot(tsb.dat[cluster != 'outlier'], aes(x = cluster, t_vaf, color = cluster))+geom_boxplot()+coord_flip()+xlab('')+ylab('')+cowplot::theme_cowplot()+
      theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())+ylim(0,1)+scale_colour_manual(values = colors)

    #Organize grid
    tsb.dat.gg = cowplot::plot_grid(tsb.dat.bar, tsb.dat.dens, nrow = 2, ncol =  1, rel_heights = c(0.25, 1))

    print(tsb.dat.gg)

    if(savePlot){
      cowplot::save_plot(filename = paste(tsb[i], '_het.pdf', sep=''), plot = tsb.dat.gg, base_height = height, base_width = width)
    }
  }
}
