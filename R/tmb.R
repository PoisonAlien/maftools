#' Estimate Tumor Mutation Burden
#' @description Estimates Tumor Mutation Burden in terms of per megabases
#' @param maf maf \code{\link{MAF}} object
#' @param captureSize capture size for input MAF in MBs. Default 50MB.
#' @param logScale Default TRUE. For plotting purpose only.
#' @return data.table with TMB for every sample
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' tmb(maf = laml)
#' @export

tmb <- function(maf, captureSize = 50, logScale = TRUE) {
  if (!is(object = maf, class2 = "MAF")) {
    stop("Input must be an MAF object")
  }

  maf.mutload <- getSampleSummary(maf)[, .(Tumor_Sample_Barcode, total)]
  maf.mutload[, total_perMB := total / captureSize]
  maf.mutload[, total_perMB_log := log10(total_perMB)]
  maf.mutload <- maf.mutload[order(total_perMB, decreasing = FALSE)]

  medload <- median(maf.mutload[, total_perMB], na.rm = TRUE)

  par(mar = c(2, 4, 2, 0))
  pointcol <- "#009688"
  medlinecol <- "#FF5722"
  if (logScale) {
    yat <- pretty(range(maf.mutload[total != 0][, total_perMB_log]))
    plot(NA, xlim = c(0, nrow(maf.mutload)), ylim = range(yat), axes = FALSE, xlab = NA, ylab = NA)
    abline(h = yat, lty = 2, lwd = 1, col = "gray90")
    abline(h = log10(medload), lty = 2, lwd = 1, col = medlinecol)
    points(x = 1:nrow(maf.mutload), y = maf.mutload$total_perMB_log, pch = 19, col = pointcol)
    title(main = "Mutation Burden", sub = paste0("Median: ", medload, "/MB"), line = 0, adj = 0)
    axis(side = 2, at = yat, las = 2)
    mtext(text = "TMB/MB (log10)", side = 2, line = 2.5)
  } else {
    yat <- pretty(range(maf.mutload$total_perMB))
    plot(NA, xlim = c(0, nrow(maf.mutload)), ylim = range(yat), axes = FALSE, xlab = NA, ylab = NA)
    abline(h = yat, lty = 2, lwd = 1, col = "gray90")
    abline(h = medload, lty = 2, lwd = 1, col = medlinecol)
    points(x = 1:nrow(maf.mutload), y = maf.mutload$total_perMB, pch = 19, col = pointcol)
    title(main = "Mutation Burden", sub = paste0("Median: ", medload, "/MB"), line = 0, adj = 0)
    axis(side = 2, at = yat, las = 2)
    mtext(text = "TMB/MB", side = 2, line = 2.5)
  }


  maf.mutload
}
