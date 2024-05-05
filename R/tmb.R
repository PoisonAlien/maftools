#' Estimate Tumor Mutation Burden
#' @description Estimates Tumor Mutation Burden in terms of per megabases
#' @param maf maf \code{\link{MAF}} object
#' @param captureRegions capture regions. Default NULL. If provided sub-sets variants within the capture regions for TMB estimation.
#' Can be a data.frame or a tsv with first three columns containing chromosome, start and end position.
#' @param captureSize capture size for input MAF in MBs. Default 50MB. Mutually exclusive with \code{captureRegions}
#' @param logScale Default TRUE. For plotting purpose only.
#' @param ignoreCNV Default TRUE. Ignores all the variants annotated as `CNV` in the `Variant_Type` column of MAF
#' @param plotType Can be "classic" or "boxplot". Set to `NA` for no plot.
#' @param pointcol Default #2c3e50
#' @return data.table with TMB for every sample
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' tmb(maf = laml)
#' @export

tmb = function(maf, captureRegions = NULL, captureSize = 50, logScale = TRUE,
               ignoreCNV = TRUE, plotType = "classic", pointcol = "#2c3e50", verbose = TRUE){

  if(!is(object = maf, class2 = "MAF")){
    stop("Input must be an MAF object")
  }

  if(!is.na(plotType)){
    plotType = match.arg(arg = plotType, choices = c("classic", "boxplot"))
  }

  if(!is.null(captureRegions)){
    if(is(object = captureRegions, class2 = "data.frame")){
      captureRegions = data.table::as.data.table(x = captureRegions)
    }else if(file.exists(captureRegions)){
      captureRegions = data.table::fread(input = captureRegions)
    }else{
      stop("captureRegions must be a tsv file or a data.frame!")
    }
    colnames(captureRegions)[1:3] = c("chrom", "start", "end")
    captureRegions[, size := end - start]
    captureSize = captureRegions[,sum(size)]/1e6
    if(verbose){
      message("Capture size:", captureSize, " Mb [from ", nrow(captureRegions), " regions]")
      message("Subsetting to provided capture regions:")
    }
    maf = maftools::subsetMaf(maf = maf, ranges = captureRegions)
  }

  if(ignoreCNV){
    if(verbose){
      message("Filtering CNV events (if any..)")
    }

    maf = subsetMaf(maf = maf, query = "!Variant_Type %in% 'CNV'", verbose = verbose)
  }

  maf.mutload = getSampleSummary(maf)[,.(Tumor_Sample_Barcode, total)]
  maf.mutload[,total_perMB := total/captureSize]
  maf.mutload[,total_perMB_log := log10(total_perMB)]
  maf.mutload = maf.mutload[order(total_perMB, decreasing = FALSE)]

  medload = median(maf.mutload[,total_perMB], na.rm = TRUE)

  if(!is.na(plotType)){
    medlinecol = "#FF5722"
    ylab = c(0.001, 0.01, 0.1, 1, 10, 100, 1000)
    yat = log10(ylab)
    par(mar = c(2, 4, 2, 0))

    if(plotType == "classic"){
      if(logScale){
        plot(NA, xlim = c(0, nrow(maf.mutload)), ylim = range(yat), axes = FALSE, xlab = NA, ylab = NA)
        #abline(h = yat, lty = 2, lwd = 1, col = 'gray90')
        grid()
        abline(h = log10(medload), lty = 2, lwd = 1, col = medlinecol)
        points(x = 1:nrow(maf.mutload), y = maf.mutload$total_perMB_log, pch = 19, col = adjustcolor(pointcol, 0.4))
        title (main = NA, sub = paste0("Median: ", round(medload, 2), "/MB"), line = 0, adj = 0)
        axis(side = 2,at = yat, las = 2, labels = ylab)
        mtext(text = "TMB/MB (log10)", side = 2, line = 2.5)
      }else{
        yat = pretty(range(maf.mutload$total_perMB))
        plot(NA, xlim = c(0, nrow(maf.mutload)), ylim = range(yat), axes = FALSE, xlab = NA, ylab = NA)
        #abline(h = yat, lty = 2, lwd = 1, col = 'gray90')
        grid()
        abline(h = medload, lty = 2, lwd = 1, col = medlinecol)
        points(x = 1:nrow(maf.mutload), y = maf.mutload$total_perMB, pch = 19, col = adjustcolor(pointcol, 0.4))
        title (main = NA, sub = paste0("Median: ", round(medload, 2), "/MB"), line = 0, adj = 0)
        axis(side = 2,at = yat, las = 2)
        mtext(text = "TMB/MB", side = 2, line = 2.5)
      }
    }else if(plotType == "boxplot"){
      if(logScale){
        ylab = c(0.01, 0.1, 1, 10, 100)
        yat = log10(ylab)
        b = boxplot(maf.mutload$total_perMB_log, border = "#2c3e50", ylim = range(yat), axes = FALSE, names = NA, ylab = "TMB (log10)", outline = FALSE, outwex = 0, col = NA)
        axis(side = 2, at = yat, labels = ylab, las = 2)
        stripchart(maf.mutload$total_perMB_log, col = adjustcolor(pointcol, 0.4), vertical = TRUE, add = TRUE, method = "jitter", pch  = 19, cex = 0.6)
        axis(side = 3, at = 1:length(b$n), labels = paste0("N: ", b$n), font =1, tick = FALSE, line = -1, cex.axis = 1)
        grid()
        title (main = NA, sub = paste0("Median: ", round(medload, 2), "/MB"), line = 0, adj = 0)
      }else{
        #ylab = c(0.01, 0.1, 1, 10, 100)
        yat = pretty(range(maf.mutload$total_perMB))
        b = boxplot(maf.mutload$total_perMB, border = "#2c3e50", ylim = range(yat), axes = FALSE, names = NA, ylab = "TMB", outline = FALSE, outwex = 0, col = NA)
        axis(side = 2, at = yat, las = 2)
        stripchart(maf.mutload$total_perMB, col = adjustcolor(pointcol, 0.4), vertical = TRUE, add = TRUE, method = "jitter", pch  = 19, cex = 0.6)
        axis(side = 3, at = 1:length(b$n), labels = paste0("N: ", b$n), font =1, tick = FALSE, line = -1, cex.axis = 1)
        grid()
        title (main = NA, sub = paste0("Median: ", round(medload, 2), "/MB"), line = 0, adj = 0)
      }
    }
  }

  maf.mutload
}

