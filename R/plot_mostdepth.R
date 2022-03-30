#'Plot results from mosdepth output
#' @param bed mosdepth output
#' @param col Colors. Default c("#95a5a6", "#7f8c8d")
#' @param sample_name sample name. Default parses from `bed`
#' @param segment Performs CBS segmentation. Default FALSE
#' @references Pedersen BS, Quinlan AR. Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics. 2018;34(5):867-868. doi:10.1093/bioinformatics/btx699
#' @export

plotMosdepth_t = function(bed = NULL, col = c("#95a5a6", "#7f8c8d"), sample_name = NULL, segment = FALSE){

  tum_cov = data.table::fread(input = bed)
  colnames(tum_cov) = c("chr", "start", "end", "doc")

  contigs = c(1:22, "X", "Y", paste0("chr", 1:22), "chrX", "chrY")
  tum_cov = tum_cov[chr %in% contigs]

  tnxy = tum_cov[chr %in% c('X', 'Y', 'chrX', 'chrY')]
  tn = tum_cov[!chr %in% c('X', 'Y', 'chrX', 'chrY')]
  tn[, chr := gsub(pattern = "chr", replacement = "", x = chr)]
  tn[, chr := as.numeric(as.character(chr))]
  tn = tn[order(chr, start)]
  med_cov = median(tn[,doc], na.rm = TRUE) #Median coverage from autosomes only

  all_depth = rbind(tn, tnxy)
  colnames(all_depth)[1:3] = c("Chromosome", "Start_Position", "End_Position")

  #Extract chromosome lengths
  chr.lens.dt = all_depth[,max(End_Position, na.rm = TRUE), .(Chromosome)]
  chr.lens = chr.lens.dt$V1
  names(chr.lens) = chr.lens.dt$Chromosome

  #chr.lens = unlist(lapply(split(all_depth, all_depth$Chromosome), function(x) max(x$End_Position, na.rm = TRUE)))

  #all_depth = .transformSegments(segmentedData = all_depth, chr.lens = chr.lens)

  if(is.null(sample_name)){
    sample_name = gsub(x = basename(bed), pattern = "\\.regions\\.bed\\.gz$", replacement = "")
  }

  #log2 median centered
  all_depth[, doc_norm := log2(doc) - log2(med_cov)]

  cn_segs = NULL
  if(segment){
    message("Running CBS segmentation:")
    #samp.name = gsub(pattern = '.denoisedCR.tsv', replacement = '', x = copynumber_file)
    cn = DNAcopy::CNA(genomdat = all_depth[,doc_norm], chrom = all_depth[,Chromosome], maploc = all_depth[,Start_Position],
                      data.type = "logratio", sampleid = sample_name, presorted = TRUE)

    cn = DNAcopy::smooth.CNA(cn)
    cn = DNAcopy::segment(cn, alpha = 0.01, nperm = 10000, p.method = 'hybrid', min.width = 5, kmax = 25, nmin = 210,
                          eta = 0.05, trim = 0.025, undo.SD = 3, undo.prune = 0.05, undo.splits = 'sdundo', verbose = 2)
    cn_segs = DNAcopy::segments.p(x = cn)
    colnames(cn_segs)[1:4] = c("Sample_Name", "Chromosome", "Start_Position", "End_Position")
    data.table::setDT(cn_segs)
    data.table::fwrite(x = cn_segs, file = paste0(sample_name, "_cbs.seg"), sep = "\t")
    message("Segments are written to: ", paste(sample_name, '_cbs.seg', sep=''))
    cn_segs = split(cn_segs, cn_segs$Chromosome)

    seg.spl.transformed = cn_segs[[1]]
    if (nrow(seg.spl.transformed) > 0) {
      seg.spl.transformed$Start_Position_updated = seg.spl.transformed$Start_Position
      seg.spl.transformed$End_Position_updated = seg.spl.transformed$End_Position
    }
    chr.lens.sumsum = cumsum(as.numeric(chr.lens))
    for (i in 2:length(cn_segs)) {
      x.seg = cn_segs[[i]]
      if (nrow(x.seg) > 0) {
        x.seg$Start_Position_updated = x.seg$Start_Position +
          chr.lens.sumsum[i - 1]
        x.seg$End_Position_updated = x.seg$End_Position +
          chr.lens.sumsum[i - 1]
      }
      seg.spl.transformed = rbind(seg.spl.transformed, x.seg,
                                  fill = TRUE)
    }
    seg.spl.transformed = split(seg.spl.transformed, seg.spl.transformed$Chromosome)[names(chr.lens)]
  }

  cols = rep(x = col, length(chr.lens))

  all_depth_spl = split(all_depth, all_depth$Chromosome)
  all_depth_spl = all_depth_spl[names(chr.lens)]

  all_depth_spl.transformed = all_depth_spl[[1]]
  if (nrow(all_depth_spl.transformed) > 0) {
    all_depth_spl.transformed$Start_Position_updated = all_depth_spl.transformed$Start_Position
    all_depth_spl.transformed$End_Position_updated = all_depth_spl.transformed$End_Position
  }
  chr.lens.sumsum = cumsum(as.numeric(chr.lens))
  for (i in 2:length(all_depth_spl)) {
    x.seg = all_depth_spl[[i]]
    if (nrow(x.seg) > 0) {
      x.seg$Start_Position_updated = x.seg$Start_Position +
        chr.lens.sumsum[i - 1]
      x.seg$End_Position_updated = x.seg$End_Position +
        chr.lens.sumsum[i - 1]
    }
    all_depth_spl.transformed = rbind(all_depth_spl.transformed, x.seg, fill = TRUE)
  }

  all_depth_spl = split(all_depth_spl.transformed, all_depth_spl.transformed$Chromosome)
  all_depth_spl = all_depth_spl[names(chr.lens)]

  message("Plotting")
  #png(filename = paste0(sample_name, ".png"), width = 1024, height = 600, bg = "white", res = 70)

  par(mfrow = c(2, 1), mar = c(3, 3, 2, 1))

  plot(NA, xlim = c(0, sum(chr.lens)), ylim = c(med_cov-50, med_cov+50), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
  temp = lapply(seq_along(all_depth_spl), function(idx){
    x = all_depth_spl[[idx]]
    points(x$Start_Position_updated, x$doc, pch = 19, cex = 0.5, col = cols[idx])
    rect(xleft = x[, Start_Position_updated][1], ybottom = med_cov,
         xright = x[,End_Position_updated][nrow(x)], ytop = med_cov)
  })
  abline(v = cumsum(as.numeric(chr.lens)), lty = 2)
  axis(side = 1, at = cumsum(as.numeric(chr.lens)), labels = names(chr.lens))
  axis(side = 2, at = pretty(c(med_cov-50, med_cov+50)), las = 2)
  title(main = "DOC")

  plot(NA, xlim = c(0, sum(chr.lens)), ylim = c(-3, 3), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
  temp = lapply(seq_along(all_depth_spl), function(idx){
    x = all_depth_spl[[idx]]
    points(x$Start_Position_updated, x$doc_norm, pch = 19, cex = 0.5, col = cols[idx])
    rect(xleft = x[, Start_Position_updated][1], ybottom = log2(med_cov),
         xright = x[,End_Position_updated][nrow(x)], ytop = log2(med_cov))

    if(segment){
      xs = seg.spl.transformed[[idx]]
      rect(xleft = xs$Start_Position_updated, ybottom = xs$seg.mean,
           xright = xs$End_Position_updated, ytop = xs$seg.mean, col = "maroon", lwd = 1, border = "maroon")
    }
  })
  abline(v = cumsum(as.numeric(chr.lens)), lty = 2)
  axis(side = 1, at = cumsum(as.numeric(chr.lens)), labels = names(chr.lens))
  axis(side = 2, at = seq(-3, 3, 1), las = 2)
  title(main = "DOC Median centered")
}
