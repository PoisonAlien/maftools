#' Compare two lollipop plots
#' @details Draws lollipop plot for a gene from two cohorts
#' @param m1 first \code{\link{MAF}} object
#' @param m2 second \code{\link{MAF}} object
#' @param gene HGNC symbol for which protein structure to be drawn.
#' @param AACol1 manually specify column name for amino acid changes in m1. Default looks for fields 'HGVSp_Short', 'AAChange' or 'Protein_Change'.
#' @param AACol2 manually specify column name for amino acid changes in m2. Default looks for fields 'HGVSp_Short', 'AAChange' or 'Protein_Change'.
#' @param m1_name name for \code{m1} cohort. optional.
#' @param m2_name name for \code{m2} cohort. optional.
#' @param m1_label Amino acid positions to label for \code{m1} cohort. If 'all', labels all variants.
#' @param m2_label Amino acid positions to label for \code{m2} cohort. If 'all', labels all variants.
#' @param refSeqID RefSeq transcript identifier for \code{gene} if known.
#' @param proteinID RefSeq protein identifier for \code{gene} if known.
#' @param labPosAngle angle for labels. Defaults to horizonal 0 degree labels. Set to 90 for vertical; 45 for diagonal labels.
#' @param labPosSize Text size for labels. Default 3
#' @param colors named vector of colors for each Variant_Classification. Default NULL.
#' @param alpha color adjustment. Default 1
#' @param axisTextSize text size for axis labels. Default 1.
#' @param pointSize size of lollipop heads. Default 1.2
#' @param domainLabelSize text size for domain labels. Default 1.
#' @param roundedRect Default FALSE. If `TRUE` domains are drawn with rounded corners. Requires \code{berryFunctions}
#' @param showDomainLabel Label domains within the plot. Default TRUE. If FALSE domains are annotated in legend.
#' @param domainBorderCol Default "black". Set to NA to remove.
#' @param legendTxtSize Default 1.
#' @param verbose Default TRUE
#' @return invisible list of domain overlaps
#' @importFrom grDevices colors colours
#' @examples
#' primary.apl <- system.file("extdata", "APL_primary.maf.gz", package = "maftools")
#' relapse.apl <- system.file("extdata", "APL_relapse.maf.gz", package = "maftools")
#' primary.apl <- read.maf(maf = primary.apl)
#' relapse.apl <- read.maf(maf = relapse.apl)
#' lollipopPlot2(m1 = primary.apl, m2 = relapse.apl, gene = "FLT3",AACol1 = "amino_acid_change", AACol2 = "amino_acid_change", m1_name = "Primary", m2_name = "Relapse")
#' @seealso \code{\link{lollipopPlot}}
#' @seealso \code{\link{mafCompare}}
#'
#'@export
#'

lollipopPlot2 = function(m1, m2, gene = NULL, AACol1 = NULL, AACol2 = NULL, m1_name = NULL, m2_name = NULL, m1_label = NULL, m2_label = NULL, refSeqID = NULL, proteinID = NULL, labPosAngle = 0, labPosSize = 0.9,
                         colors = NULL, alpha = 1, axisTextSize = c(1, 1), pointSize = 1.2, roundedRect = TRUE, showDomainLabel = TRUE, domainBorderCol = "black", domainLabelSize = 1, legendTxtSize = 1, verbose = TRUE){

  if(is.null(gene)){
    stop('Please provide a gene name.')
  }else{
    if(verbose){
      message(paste0("Gene: ", gene))
    }
  }

  m1.lp = get_lp_data(maf = m1, geneID = gene, AACol = AACol1, refSeqID = refSeqID, proteinID = proteinID, verbose = verbose)
  m2.lp = get_lp_data(maf = m2, geneID = gene, AACol = AACol2, refSeqID = refSeqID, proteinID = proteinID, verbose = verbose)

  if(is.null(m1.lp) & is.null(m2.lp)){
    stop("At least one of the maf should contain mutations")
  }

  if(is.null(colors)){
    col = get_vcColors()
  }else{
    col = colors
  }

  if(!is.null(m1.lp) & !is.null(m2.lp)){

    m1.y.at = m1.lp[[3]]
    m1.y.at.lbl = m1.lp[[4]]
    m1.mut.per = m1.lp[[5]]

    m2.y.at = m2.lp[[3]]
    m2.y.at.lbl = m2.lp[[4]]
    m2.mut.per = m2.lp[[5]]

    prot = m2.lp[[6]]
    prot.len = m2.lp[[2]]
  }

  if(is.null(m1.lp)){
    m1.y.at = c(1, 5)
    m1.y.at.lbl = m1.y.at
    m1.mut.per = "0%"

    prot = m2.lp[[6]]
    prot.len = m2.lp[[2]]
    m2.y.at = m2.lp[[3]]
    m2.y.at.lbl = m2.lp[[4]]
    m2.mut.per = m2.lp[[5]]
  }

  if(is.null(m2.lp)){
    m2.y.at = c(1, 5)
    m2.y.at.lbl = m2.y.at
    m2.mut.per = "0%"

    prot = m1.lp[[6]]
    prot.len = m1.lp[[2]]
    m1.y.at = m1.lp[[3]]
    m1.y.at.lbl = m1.lp[[4]]
    m1.mut.per = m1.lp[[5]]
  }

  domains = unique(prot[,Label])
  domain_cols = get_domain_cols()

  if(length(domains) > length(domain_cols)){
    domain_cols = sample(colours(), size = length(domains), replace = FALSE)
  }

  domain_cols = domain_cols[1:length(domains)]
  #domain_cols = grDevices::adjustcolor(col = domain_cols, alpha.f = 0.75)
  names(domain_cols) = domains

  col1 = col2 = NULL
  if(!is.null(m1.lp)){
    col1 = col[as.character(m1.lp[[1]][,Variant_Classification])]
    names_temp = names(col1)
    col1 = grDevices::adjustcolor(col = col1, alpha.f = alpha)
    names(col1) = names_temp
    rm(names_temp)
  }

  if(!is.null(m2.lp)){
    col2 = col[as.character(m2.lp[[1]][,Variant_Classification])]
    names_temp = names(col2)
    col2 = grDevices::adjustcolor(col = col2, alpha.f = alpha)
    names(col2) = names_temp
    rm(names_temp)
  }

  lo = matrix(data = c(1, 1, 2, 2), nrow = 2, byrow = TRUE)
  graphics::layout(mat = lo, heights = c(5, 2))

  par(mar = c(2, 2.5, 2, 1))
  plot(0, 0, xlim = c(0, max(prot.len, na.rm = TRUE)+10), ylim = c(-6.5, 6.5), axes = FALSE, pch = NA, xlab = "", ylab = "")
  rect(xleft = 0, ybottom = -0.8, xright = max(prot.len), ytop = 0.8, col = "#95a5a6", border = domainBorderCol)
  # axis(side = 1, at = m1.lp[[2]], labels = m1.lp[[2]], lwd = 1.2, font = 1,
  #      cex.axis = axisTextSize[1], line = -0.4)
  axis(side = 2, at = m1.y.at, labels = m1.y.at.lbl, lwd = 1.2, font = 1, las = 2,
       cex.axis = axisTextSize[2])
  axis(side = 2, at = -m2.y.at, labels = m2.y.at.lbl, lwd = 1.2, font = 1, las = 2,
       cex.axis = axisTextSize[2])

  if(!is.null(m1.lp)){
    segments(x0 = m1.lp[[1]][,pos2], y0 = 0.8, x1 = m1.lp[[1]][,pos2], y1 = m1.lp[[1]][,count2-0.03], lwd = 1.2, col = "gray70")
    points(x = m1.lp[[1]][,pos2], y = m1.lp[[1]][,count2], col = col1, pch = 16, cex = pointSize)
  }

  if(!is.null(m2.lp)){
    segments(x0 = m2.lp[[1]][,pos2], y0 = -0.8, x1 = m2.lp[[1]][,pos2], y1 = -m2.lp[[1]][,count2-0.03], lwd = 1.2, col = "gray70")
    points(x = m2.lp[[1]][,pos2], y = -m2.lp[[1]][,count2], col = col2, pch = 16, cex = pointSize)
  }

  col_uniq = unique(c(names(col1), names(col2)))
  col = c(col1, col2)[col_uniq]

  prot[, domainCol := domain_cols[prot[, Label]]]
  if(roundedRect){
    if(requireNamespace("berryFunctions", quietly = TRUE)){
      for(i in 1:nrow(prot)){
        berryFunctions::roundedRect(xleft = prot[i,Start], ybottom = -1, xright = prot[i,End], ytop = 1, col = prot[i, domainCol], border = domainBorderCol, rounding = 0.08)
      }
    }else{
      #warning("Package berryFunctions needed for roundedRect to work. Please install it and try again.")
      rect(xleft = prot[,Start], ybottom = -1, xright = prot[,End], ytop = 1, col = prot[,domainCol], border = domainBorderCol)
    }
  }else{
    rect(xleft = prot[,Start], ybottom = -1, xright = prot[,End], ytop = 1, col = prot[,domainCol], border = domainBorderCol)
  }
  #rect(xleft = prot[,Start], ybottom = -1, xright = prot[,End], ytop = 1, col = prot[,domainCol], border = NA)

  prot$pos = rowMeans(x = prot[,.(Start, End)])
  if(showDomainLabel){
    prot = prot[!duplicated(Label)]
    text(y = 0, x = prot$pos, labels = prot$Label, font = 3, cex = domainLabelSize)
  }


  text(x = max(prot.len)+12, y = 0.5, labels = paste0(max(prot.len), " aa"), font = 1, adj = 1, srt = 90)


  if(!is.null(m1_name)){
    mtext(text = paste0(m1_name, " [", m1.mut.per, "; N = ", m1@summary[3, summary], "]"), side = 3, line = 0.2, adj = 0, font = 4)
  }else{
    mtext(text = paste0(m1.mut.per, "; N = ", m1@summary[3, summary]), side = 3, line = 0.2, adj = 0, font = 3)
  }

  mtext(text = paste0(gene, ": ", unique(prot[,refseq.ID])), side = 3, line = 0.2, adj = 1, font = 1)

  if(!is.null(m2_name)){
    mtext(text = paste0(m2_name, " [", m2.mut.per, "; N = ", m2@summary[3, summary], "]"), side = 1, line = 0.2, adj = 0, font = 4)
  }else{
    mtext(text = paste0(m2.mut.per, "; N = ", m2@summary[3, summary]), side = 1, line = 0.2, adj = 0, font = 3)
  }

  if(!is.null(m1_label)){
    if(!is.null(m1.lp)){
      m1_label = label_pos(prot.snp.sumamry = m1.lp[[1]], labelPos = m1_label)
      if(!is.null(m1_label)){
        text(x = m1_label[,pos2], y = m1_label[,count2+0.45], labels = m1_label[,conv],
             font = 1, srt = labPosAngle, cex = labPosSize)
      }
    }
  }

  if(!is.null(m2_label)){
    if(!is.null(m2.lp)){
      m2_label = label_pos(prot.snp.sumamry = m2.lp[[1]], labelPos = m2_label)
      if(!is.null(m2_label)){
        text(x = m2_label[,pos2], y = -m2_label[,count2+0.60], labels = m2_label[,conv],
             font = 1, srt = labPosAngle, cex = labPosSize)
      }
    }
  }

  par(mar = c(5, 0.5, 0, 0), xpd = TRUE)

  plot(NULL,ylab='',xlab='', xlim=0:1, ylim=0:1, axes = FALSE)
  lep = legend("topleft", legend = names(col), col = col,  bty = "n", border=NA,
                 xpd = TRUE, text.font = 1, pch = 16, xjust = 0, yjust = 0,
                 cex = legendTxtSize, y.intersp = 1.5, x.intersp = 1,
                 pt.cex = 1.2 * legendTxtSize, ncol = ceiling(length(col)/4))

    x_axp = 0+lep$rect$w

    if(!showDomainLabel){
      if(length(domain_cols) <= 4){
        n_col = 1
      }else{
        n_col = (length(domain_cols) %/% 4)+1
      }

      lep = legend(x = x_axp, y = 1, legend = names(domain_cols),
                   col = domain_cols, border = NA,
                   ncol= n_col, pch = 15, xpd = TRUE, xjust = 0, bty = "n",
                   cex = legendTxtSize, title = "Domains",
                   title.adj = 0, pt.cex = 1.2 * legendTxtSize)
    }

    invisible(list(M1 = m1.lp[[7]], M2 = m2.lp[[7]]))

}
