#' Plot oncogenic pathways
#' @references Sanchez-Vega F, Mina M, Armenia J, Chatila WK, Luna A, La KC, Dimitriadoy S, Liu DL, Kantheti HS, Saghafinia S et al. 2018. Oncogenic Signaling Pathways in The Cancer Genome Atlas. Cell 173: 321-337 e310
#'
#' @details
#' Draws pathway burden123
#' @param maf an \code{\link{MAF}} object
#' @param pathlist Output from  \code{\link{pathways}}
#' @param pathnames Names of the pathways to be drawn. Default NULL, plots everything from input `pathlist`
#' @param removeNonMutated Default FALSE
#' @param fontSize Default 1
#' @param showTumorSampleBarcodes logical to include sample names.
#' @param SampleNamefontSize font size for sample names. Default 0.6
#' @param sampleOrder Manually speify sample names for oncolplot ordering. Default NULL.
#' @param mar margins Default c(4, 6, 2, 3). Margins to bottom, left, top and right respectively
#' @seealso \code{\link{pathways}}
#' @export
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' p <- pathways(maf = laml)
#' plotPathways(maf = laml, pathlist = p)
plotPathways = function(maf = NULL, pathlist = NULL, pathnames = NULL,
                                 removeNonMutated = FALSE, fontSize = 1, showTumorSampleBarcodes = FALSE, sampleOrder = NULL,
                        SampleNamefontSize = 0.6, mar = c(4, 6, 2, 3)){

  if(any(is.null(pathlist), is.null(maf))){
    stop("Misssing required input maf or pathlist")
  }

  pws = attr(pathlist, "genes")
  if(!is.null(pathnames)){
    pathnames = match.arg(arg = pathnames, choices = names(pws), several.ok = TRUE)
  }
  #print(pws)

  totSamps = as.numeric(maf@summary[3,summary])
  tsbs = levels(getSampleSummary(x = maf)[,Tumor_Sample_Barcode])

  genes = unique(unlist(pws, use.names = FALSE))
  path_mat = createOncoMatrix(m = maf, g = genes)

  nm = path_mat$numericMatrix
  genes.missing = genes[!genes %in% rownames(nm)]
  if(length(genes.missing) > 0){
    genes.missing.numat = matrix(data = 0, ncol = ncol(nm), nrow = length(genes.missing))
    rownames(genes.missing.numat) = genes.missing
    nm = rbind(nm, genes.missing.numat)
  }

  genes = rownames(nm)

  if(!removeNonMutated){
    tsb.include = matrix(data = 0, nrow = length(genes), ncol = length(tsbs[!tsbs %in% colnames(nm)]))
    colnames(tsb.include) = tsbs[!tsbs %in% colnames(nm)]
    rownames(tsb.include) = rownames(nm)
    nm = cbind(nm, tsb.include)
  }

  nm = t(apply(nm, 2, rev))
  nm[nm == 0] = NA
  nm[!is.na(nm)] = 1

  collapse = TRUE
  if(collapse){
    pws_mat = lapply(pws, function(p){
      nm_p = nm[,p, drop = FALSE]
      nm_p[is.na(nm_p)] = 0
      nm_p = rowSums(nm_p, na.rm = TRUE)
      as.matrix(x = nm_p)
    })
    pws_mat = do.call(what = "cbind", pws_mat)
    colnames(pws_mat) = names(pws)
    pws_mat[pws_mat > 0] = 1

    if(!is.null(pathnames)){
      pws_mat = pws_mat[,pathnames, drop = FALSE]
    }

    if(ncol(pws_mat) == 1){
      pws_mat = pws_mat[,names(sort(colSums(pws_mat), decreasing = TRUE)), drop = FALSE]
      pws_mat = pws_mat[order(pws_mat[,1], decreasing = TRUE),, drop = FALSE]
      pws_mat[pws_mat == 0] = NA
    }else{
      pws_mat = pws_mat[,names(sort(colSums(pws_mat), decreasing = TRUE)), drop = FALSE]
      pws_mat = t(pws_mat[do.call(order, c(as.list(as.data.frame(pws_mat)), decreasing = TRUE)), ])
      pws_mat[pws_mat == 0] = NA
      pws_mat = t(apply(pws_mat, 2, rev))
    }

    pws_mat_bg = pws_mat
    pws_mat_bg[is.na(pws_mat_bg)] = 1
  }

  if(!is.null(sampleOrder)){
    sampleOrder = as.character(sampleOrder)
    sampleOrder = intersect(sampleOrder, rownames(pws_mat))

    if(length(sampleOrder) == 0){
      stop("None of the provided samples are present in the input MAF")
    }
    pws_mat = pws_mat[sampleOrder, ,drop = FALSE]
  }

  pw_pct = paste0(round((colSums(pws_mat, na.rm = TRUE)/totSamps) * 100, digits = 2), "%")
  par(mar = mar)

  image(x = 1:nrow(pws_mat_bg), y = 1:ncol(pws_mat_bg), z = pws_mat_bg, axes = FALSE, xaxt="n", yaxt="n", xlab="", ylab="", col = "#ecf0f1")
  image(x = 1:nrow(pws_mat), y = 1:ncol(pws_mat), z = pws_mat, axes = FALSE, xaxt="n", yaxt="n", xlab="", ylab="", col = "#34495e", add = TRUE) #col = "#FC8D62"
  abline(h = (1:ncol(pws_mat)) + 0.5, col = "white")
  abline(v = (1:nrow(pws_mat)) + 0.5, col = "white")
  mtext(text = colnames(pws_mat), side = 2, at = 1:ncol(pws_mat), line = 0.4, cex = fontSize, las = 2)
  mtext(text = pw_pct, side = 4, at = 1:ncol(pws_mat), line = 0.4, cex = fontSize, las = 2)

  if(showTumorSampleBarcodes){
    text(y = rep(0, nrow(pws_mat)), x = 1:nrow(pws_mat),
         labels = rownames(pws_mat), srt = 45, font = 1, cex = SampleNamefontSize, adj = 1, xpd = TRUE)
  }

  invisible(rownames(pws_mat))
}
