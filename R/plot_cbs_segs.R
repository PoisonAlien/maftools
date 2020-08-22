#' Plots segmented copy number data.
#' @details this function takes segmented copy number data and plots it. If MAF object is specified, all mutations are highlighted on the plot.
#' @param cbsFile CBS segmented copy number file. Column names should be Sample, Chromosome, Start, End, Num_Probes and Segment_Mean (log2 scale).
#' @param maf optional \code{\link{MAF}}
#' @param tsb If segmentation file contains many samples (as in gistic input), specify sample name here.
#' Defualt plots head 1 sample. Set 'ALL' for plotting all samples. If you are maping maf, make sure sample names in
#' Sample column of segmentation file matches to those Tumor_Sample_Barcodes in MAF.
#' @param savePlot If true plot is saved as pdf.
#' @param ylims Default NULL
#' @param seg_size Default 0.1
#' @param width width of plot
#' @param height height of plot
#' @param genes If given and maf object is specified, maps all mutataions from maf onto segments. Default NULL
#' @param ref.build Reference build for chromosome sizes. Can be hg18, hg19 or hg38. Default hg19.
#' @param writeTable If true and if maf object is specified, writes plot data with each variant and its corresponding copynumber to an output file.
#' @param removeXY don not plot sex chromosomes.
#' @param color Manually specify color scheme for chromosomes. Default NULL. i.e, aletrnating Gray70 and midnightblue
#' @return Draws plot
#' @export
#' @examples
#' tcga.ab.009.seg <- system.file("extdata", "TCGA.AB.3009.hg19.seg.txt", package = "maftools")
#' plotCBSsegments(cbsFile = tcga.ab.009.seg)
plotCBSsegments <- function(cbsFile = NULL, maf = NULL, tsb = NULL, savePlot = FALSE, ylims = NULL, seg_size = 0.1,
                            width = 6, height = 3, genes = NULL, ref.build = "hg19", writeTable = FALSE, removeXY = FALSE, color = NULL) {
  if (is.null(cbsFile)) {
    stop("Missing segmentation file!")
  }

  # Read segmentation file and change chromosome names
  seg <- readSegs(seg = cbsFile)

  if (removeXY) {
    seg <- seg[!Chromosome %in% c("23", "24")]
  }

  seg <- seg[order(as.numeric(Chromosome))]
  data.table::setkey(x = seg, Chromosome, Start_Position, End_Position)

  # If user doesn't provide sample name
  if (is.null(tsb)) {
    # Use top 1 sample for simplicity
    tsb <- unique(as.character(seg[, Sample]))[1]
    message("No 'tsb' specified, plot head 1 sample. Set tsb='ALL' to plot all samples.")
  } else if (tsb == "ALL") {
    # Number of unique samples in segmentation file
    tsb <- unique(as.character(seg[, Sample]))
  }

  for (i in 1:length(tsb)) {
    if (savePlot) {
      pdf(file = paste(tsb[i], "_segPlot.pdf", sep = ""), width = width, height = height, bg = "white", paper = "special")
    }

    plotCBS(
      segData = seg, tsb = tsb[i], build = ref.build, chr.colors = color,
      y_lims = ylims, rect_size = seg_size
    )

    if (!is.null(maf)) {
      tsb.mapped <- mapMutsToSegs(seg = seg, maf = maf, tsb = tsb[i], build = ref.build)
      if (writeTable) {
        write.table(
          x = tsb.mapped[, .(Hugo_Symbol, Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode, Segment_Start, Segment_End, Segment_Mean, CN)],
          file = paste(tsb[i], "_segData.tsv", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE
        )
      }
    }

    if (!is.null(genes)) {
      if (is.null(maf)) {
        warning("Missing maf file. Skipping gene mapping..", immediate. = TRUE)
      } else {
        # Map mutations to segments
        tsb.mapped.cn <- tsb.mapped[Hugo_Symbol %in% genes]

        if (nrow(tsb.mapped.cn) > 0) {
          text(
            x = tsb.mapped.cn$Start_Position_updated,
            y = ifelse(test = tsb.mapped.cn$Segment_Mean > 0, yes = tsb.mapped.cn$Segment_Mean + 1, no = tsb.mapped.cn$Segment_Mean - 1),
            labels = as.character(tsb.mapped.cn$Hugo_Symbol), cex = 0.8
          )
          segments(
            x0 = tsb.mapped.cn$Start_Position_updated, y0 = tsb.mapped.cn$Segment_Mean,
            x1 = tsb.mapped.cn$Start_Position_updated,
            y1 = ifelse(test = tsb.mapped.cn$Segment_Mean > 0, yes = tsb.mapped.cn$Segment_Mean + 0.9, no = tsb.mapped.cn$Segment_Mean - 0.9),
            lty = 1
          )
        } else {
          message(paste0("No mutations observed from given gene list for sample ", tsb[i]))
        }
      }
    }
    if (savePlot) {
      dev.off()
    }
  }
  return(invisible(NULL))
}
