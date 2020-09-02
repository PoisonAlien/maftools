#' Draw two oncoplots side by side for cohort comparision.
#' @details Draws two oncoplots side by side to display difference between two cohorts.
#'
#' @param m1 first \code{\link{MAF}} object
#' @param m2 second \code{\link{MAF}} object
#' @param genes draw these genes. Default plots top 5 mutated genes from two cohorts.
#' @param clinicalFeatures1 columns names from `clinical.data` slot of m1 \code{MAF} to be drawn in the plot. Dafault NULL.
#' @param clinicalFeatures2 columns names from `clinical.data` slot of m2 \code{MAF} to be drawn in the plot. Dafault NULL.
#' @param annotationColor1 list of colors to use for `clinicalFeatures1` Default NULL.
#' @param annotationColor2 list of colors to use for `clinicalFeatures2` Default NULL.
#' @param sortByAnnotation1 logical sort oncomatrix (samples) by provided `clinicalFeatures1`. Sorts based on first `clinicalFeatures1`.  Defaults to FALSE. column-sort
#' @param sortByAnnotation2 same as above but for m2
#' @param sampleOrder1 Manually speify sample names in m1 for oncolplot ordering. Default NULL.
#' @param sampleOrder2 Manually speify sample names in m2 for oncolplot ordering. Default NULL.
#' @param additionalFeature1 a vector of length two indicating column name in the MAF and the factor level to be highlighted.
#' @param additionalFeaturePch1 Default 20
#' @param additionalFeatureCol1 Default "white"
#' @param additionalFeatureCex1 Default 0.9
#' @param additionalFeature2 a vector of length two indicating column name in the MAF and the factor level to be highlighted.
#' @param additionalFeaturePch2 Default 20
#' @param additionalFeatureCol2 Default "white"
#' @param additionalFeatureCex2 Default 0.9
#' @param sepwd_genes1 Default 0.5
#' @param sepwd_samples1 Default 0.5
#' @param sepwd_genes2 Default 0.5
#' @param sepwd_samples2 Default 0.5
#' @param annotationFontSize font size for annotations Default 1.2
#' @param colors named vector of colors for each Variant_Classification.
#' @param removeNonMutated Logical. If \code{TRUE} removes samples with no mutations in the oncoplot for better visualization. Default \code{TRUE}.
#' @param m1Name optional name for first cohort
#' @param m2Name optional name for second cohort
#' @param geneNamefont font size for gene names. Default 1
#' @param showSampleNames whether to show sample names. Defult FALSE.
#' @param barcode_mar Margin width for sample names. Default 1
#' @param gene_mar Margin width for gene names. Default 3
#' @param SampleNamefont font size for sample names. Default 0.5
#' @param anno_height Height of clinical margin. Default 2
#' @param legend_height Height of legend margin. Default 4
#' @param legendFontSize font size for legend. Default 1.2
#' @param titleFontSize font size for title. Default 1.5
#' @param keepGeneOrder force the resulting plot to use the order of the genes as specified. Default FALSE
#' @param bgCol Background grid color for wild-type (not-mutated) samples. Default gray - "#CCCCCC"
#' @param borderCol border grid color for wild-type (not-mutated) samples. Default 'white'
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
#'
coOncoplot = function(m1, m2, genes = NULL, m1Name = NULL, m2Name = NULL,
                       clinicalFeatures1 = NULL, clinicalFeatures2 = NULL,
                       annotationColor1 = NULL, annotationColor2 = NULL, annotationFontSize = 1.2,
                       sortByAnnotation1 = FALSE, sortByAnnotation2 = FALSE,
                      sampleOrder1 = NULL, sampleOrder2 = NULL,
                      additionalFeature1 = NULL, additionalFeaturePch1 = 20, additionalFeatureCol1 = "white", additionalFeatureCex1 = 0.9,
                      additionalFeature2 = NULL, additionalFeaturePch2 = 20, additionalFeatureCol2 = "white", additionalFeatureCex2 = 0.9,
                      sepwd_genes1 = 0.5, sepwd_samples1 = 0.5, sepwd_genes2 = 0.5, sepwd_samples2 = 0.5,
                       colors = NULL, removeNonMutated = TRUE, anno_height = 2, legend_height = 4,
                       geneNamefont = 0.8, showSampleNames = FALSE, SampleNamefont = 0.5, barcode_mar = 1,gene_mar = 3,
                       legendFontSize = 1.2, titleFontSize = 1.5, keepGeneOrder=FALSE,
                       bgCol = "#CCCCCC", borderCol = "white"){

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

  m1Name = paste(m1Name, ' (N = ' , m1.sampleSize, ')',sep = '')

  if(is.null(m2Name)){
    m2Name = 'M2'
  }

  m2Name = paste(m2Name, ' (N = ' , m2.sampleSize, ')',sep = '')

  if(is.null(colors)){
    vc_col = get_vcColors(websafe = FALSE)
  }else{
    vc_col = colors
  }

  m12_annotation_colors = NULL
  if(!is.null(clinicalFeatures1) & !is.null(clinicalFeatures2)){
    if(is.null(annotationColor1) & is.null(annotationColor2)){
      m12_annotation_colors = get_m12_annotation_colors(a1 = m1, a1_cf = clinicalFeatures1,
                                                        a2 = m2, a2_cf = clinicalFeatures2)
      annotationColor1 = m12_annotation_colors
      annotationColor2 = m12_annotation_colors
    }
  }

  #Get matrix dimensions and legends to adjust plot graphics::layout
  nm1 = print_mat(maf = m1, genes = genes, removeNonMutated = removeNonMutated,
                  test = TRUE, colors = colors)
  #print(nm1)
  nm1_ncol = ncol(nm1[[1]])
  nm1_vc_cols = nm1[[2]]

  nm2 = print_mat(maf = m2, genes = genes, removeNonMutated = removeNonMutated,
                  test = TRUE, colors = colors)
  nm2_ncol = ncol(nm2[[1]])
  nm2_vc_cols = nm2[[2]]

  if(!is.null(clinicalFeatures1) || !is.null(clinicalFeatures2)){
    mat_lo = matrix(data = c(1,3,5,2,4,6,7,7,7), nrow = 3, ncol = 3, byrow = TRUE)
    mat_lo = graphics::layout(mat = mat_lo,
                    widths = c(6 * (nm1_ncol/nm2_ncol), 1.5, 6), heights = c(12, anno_height, legend_height))
  }else{
    mat_lo = matrix(data = c(1,2,3,4,4,4), nrow = 2, ncol = 3, byrow = TRUE)
    mat_lo = graphics::layout(mat = mat_lo,
                    widths = c(6 * (nm1_ncol/nm2_ncol), 1, 6), heights = c(12, legend_height))
  }

  #Plot first oncoplot
  m1_legend = print_mat(maf = m1, genes = genes, removeNonMutated = removeNonMutated,
                        clinicalFeatures = clinicalFeatures1, colors = colors,
                        annotationColor = annotationColor1, barcode_size = SampleNamefont,
                        sortByAnnotation = sortByAnnotation1, fontSize = geneNamefont,
                        title = m1Name, title_size = titleFontSize,
                        showBarcodes = showSampleNames, bgCol = bgCol, borderCol = borderCol,
                        additionalFeature = additionalFeature1, additionalFeaturePch = additionalFeaturePch1,
                        additionalFeatureCex = additionalFeatureCex1, additionalFeatureCol = additionalFeatureCol1,
                        sepwd_genes = sepwd_genes1, sepwd_samples = sepwd_samples1, barcodemar = barcode_mar, genemar = gene_mar, sampleOrder = sampleOrder1)


  if(is.null(clinicalFeatures1) & !is.null(clinicalFeatures2)){
    plot.new()
  }

  if(showSampleNames){
    par(mar = c(barcode_mar, 0, 3, 0))
  }else{
    par(mar = c(1, 0, 3, 0))
  }

  #Plot gene names
  nm = matrix(data = 1, nrow = 1, ncol = length(genes))
  image(x = 1:nrow(nm), y = 1:ncol(nm), z = nm, axes = FALSE, xaxt="n", yaxt="n",
        xlab="", ylab="", col = "white")
  mtext(text = rev(genes), side = 2, adj = 0.5, at = 1:ncol(nm),
        font = 3, line = -2, cex = geneNamefont, las = 2)

  if(!is.null(clinicalFeatures1) || !is.null(clinicalFeatures2)){
    plot.new()
  }

  m2_legend = print_mat(maf = m2, genes = genes, removeNonMutated = removeNonMutated,
                        clinicalFeatures = clinicalFeatures2, colors = colors,
                        annotationColor = annotationColor2, barcode_size = SampleNamefont,
                        sortByAnnotation = sortByAnnotation2, fontSize = geneNamefont,
                        title = m2Name, title_size = titleFontSize, plot2 = TRUE,
                        showBarcodes = showSampleNames, bgCol = bgCol, borderCol = borderCol,
                        additionalFeature = additionalFeature2, additionalFeaturePch = additionalFeaturePch2,
                        additionalFeatureCex = additionalFeatureCex2, additionalFeatureCol = additionalFeatureCol2,
                        sepwd_genes = sepwd_genes2, sepwd_samples = sepwd_samples2, barcodemar = barcode_mar, genemar = gene_mar, sampleOrder = sampleOrder2)

  if(!is.null(clinicalFeatures1) & is.null(clinicalFeatures2)){
    plot.new()
  }

  vc_legend = list(nm1_vc_cols, nm2_vc_cols)
  vc_legend = unlist(vc_legend)[unique(names(unlist(vc_legend)))]
  #vc_legend = vc_col[vc_legend]
  vc_legend = vc_legend[!is.na(vc_legend)]

  if(is.null(m12_annotation_colors)){
    anno_legend = c(m1_legend, m2_legend)
  }else{
    anno_legend = m12_annotation_colors
  }

  par(mar = c(1, 1, 0, 0), xpd = TRUE)

  vc_pch = rep(15, length(vc_legend))
  if(!is.null(additionalFeature1)){
    vc_legend = c(vc_legend, "gray70")
    names(vc_legend)[length(vc_legend)] = paste(additionalFeature1, collapse = ":")
    vc_pch = c(vc_pch, additionalFeaturePch1)
  }

  if(!is.null(additionalFeature2)){
    vc_legend = c(vc_legend, "gray70")
    names(vc_legend)[length(vc_legend)] = paste(additionalFeature2, collapse = ":")
    vc_pch = c(vc_pch, additionalFeaturePch2)
  }

  plot(NULL,ylab='',xlab='', xlim=0:1, ylim=0:1, axes = FALSE)
  lep = legend("topleft", legend = names(vc_legend),
               col = vc_legend, border = NA, bty = "n",
               ncol= 2, pch = vc_pch, xpd = TRUE, xjust = 0, yjust = 0, cex = legendFontSize)

  x_axp = 0+lep$rect$w

  if(!is.null(anno_legend)){

    for(i in 1:length(anno_legend)){
      #x = unique(annotation[,i])
      x = anno_legend[[i]]

      if(length(x) <= 4){
        n_col = 1
      }else{
        n_col = (length(x) %/% 4)+1
      }

      lep = legend(x = x_axp, y = 1, legend = names(x),
                   col = x, border = NA,
                   ncol= n_col, pch = 15, xpd = TRUE, xjust = 0, bty = "n",
                   cex = annotationFontSize, title = names(anno_legend)[i], title.adj = 0)
      x_axp = x_axp + lep$rect$w

    }
  }
  #title()
}
