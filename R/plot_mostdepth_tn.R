#'Plot results from mosdepth output for Tumor/Normal pair
#' @param t_bed mosdepth output from tumor
#' @param n_bed mosdepth output from matched normal
#' @param segment Performs CBS segmentation. Default TRUE
#' @param sample_name sample name. Default parses from `t_bed`
#' @param col Colors. Default c("#95a5a6", "#7f8c8d")
#' @references Pedersen BS, Quinlan AR. Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics. 2018;34(5):867-868. doi:10.1093/bioinformatics/btx699
#' @export

plotMosdepth = function(t_bed = NULL, n_bed = NULL, segment = TRUE, sample_name = NULL, col = c("#95a5a6", "#7f8c8d")){

  contigs = c(1:22, "X", "Y", paste0("chr", 1:22), "chrX", "chrY")

  if(is.null(sample_name)){
    sample_name = gsub(x = basename(t_bed), pattern = "\\.regions\\.bed\\.gz$", replacement = "")
  }

  # if(is.null(plot_file)){
  #   plot_file = sample_name
  # }

  dat = lapply(X = c(t_bed, n_bed), function(x){
    x = data.table::fread(input = x)
    colnames(x) = c("chr", "start", "end", "doc")
    x[, chr := gsub(pattern = "chr", replacement = "", x = chr)]
    colnames(x)[1:3] = c("Chromosome", "Start_Position", "End_Position")
    x = x[Chromosome %in% contigs]
    x
  })


  names(dat) = c("tumor", "normal")
  dat = merge(dat$tumor, dat$normal, by = c("Chromosome", "Start_Position", 'End_Position'), suffixes = c("_t", "_n"))

  dat_xy = dat[Chromosome %in% c('X', 'Y', 'chrX', 'chrY')]
  datn = dat[!Chromosome %in% c('X', 'Y', 'chrX', 'chrY')]
  datn[, Chromosome := as.numeric(as.character(Chromosome))]
  datn = datn[order(Chromosome, Start_Position)]
  dat = rbind(datn, dat_xy)

  #Get chr lengths
  chr.lens.dt = dat[,max(End_Position, na.rm = TRUE), .(Chromosome)]
  chr.lens = chr.lens.dt$V1
  names(chr.lens) = chr.lens.dt$Chromosome

  map_ratio = sum(dat$doc_t, na.rm = TRUE)/sum(dat$doc_n, na.rm = TRUE)
  message("Coverage ratio T/N: ", round(map_ratio, digits = 3))
  dat$doc_n = dat$doc_n * map_ratio

  dat[, logR := log2(doc_t+1) - log2(doc_n+1)]

  cols = rep(x = col, length(chr.lens))

  all_depth_spl = split(dat, dat$Chromosome)[names(chr.lens)]

  seg.spl.transformed = all_depth_spl[[1]]
  if (nrow(seg.spl.transformed) > 0) {
    seg.spl.transformed$Start_Position_updated = seg.spl.transformed$Start_Position
    seg.spl.transformed$End_Position_updated = seg.spl.transformed$End_Position
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
    seg.spl.transformed = rbind(seg.spl.transformed, x.seg,
                                fill = TRUE)
  }

  all_depth_spl = split(seg.spl.transformed, seg.spl.transformed$Chromosome)[names(chr.lens)]

  cn_segs = NULL
  if(segment){
    message("Running CBS segmentation:")
    #samp.name = gsub(pattern = '.denoisedCR.tsv', replacement = '', x = copynumber_file)
    cn = DNAcopy::CNA(genomdat = data.table::rbindlist(l = all_depth_spl)[,logR], chrom = data.table::rbindlist(l = all_depth_spl)[,Chromosome], maploc = data.table::rbindlist(l = all_depth_spl)[,Start_Position],
                      data.type = "logratio", sampleid = sample_name, presorted = TRUE)

    cn = DNAcopy::smooth.CNA(cn)
    cn = DNAcopy::segment(cn, alpha = 0.01, nperm = 10000, p.method = 'hybrid', min.width = 5, kmax = 25, nmin = 210,
                          eta = 0.05, trim = 0.025, undo.SD = 3, undo.prune = 0.05, undo.splits = 'sdundo', verbose = 2)
    cn_segs = DNAcopy::segments.p(x = cn)
    colnames(cn_segs)[1:4] = c("Sample_Name", "Chromosome", "Start_Position", "End_Position")
    data.table::setDT(cn_segs)
    data.table::fwrite(x = cn_segs, file = paste0(sample_name, "_cbs.seg"), sep = "\t")
    message("Segments are written to: ", paste(sample_name, '_cbs.seg', sep=''))
    cn_segs = split(cn_segs, cn_segs$Chromosome)[names(all_depth_spl)]

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

  message("Plotting")
  #png(filename = paste0(sample_name, ".png"), width = 1024, height = 600, bg = "white")
  par(mar = c(4, 4, 3, 1))
  plot(NA, xlim = c(0, sum(chr.lens)), ylim = c(-3, 3), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
  temp = lapply(seq_along(all_depth_spl), function(idx){
    message("  Chromosome: ", names(all_depth_spl)[idx])
    x = all_depth_spl[[idx]]
    points(x$Start_Position_updated, x$logR, pch = 19, cex = 0.5, col = cols[idx])
    if(segment){
      xs = seg.spl.transformed[[idx]]
      rect(xleft = xs$Start_Position_updated, ybottom = xs$seg.mean,
          xright = xs$End_Position_updated, ytop = xs$seg.mean, col = "maroon", lwd = 1, border = "maroon")
    }
  })
  abline(v = cumsum(as.numeric(chr.lens)), lty = 2, col = "gray70")
  axis(side = 1, at = cumsum(as.numeric(chr.lens)), labels = names(chr.lens))
  axis(side = 2, at = seq(-3, 3, 1), las = 2)
  mtext(text = "logR", side = 2, line = 2)
  title(main = sample_name, adj = 0)
  #dev.off()
}
