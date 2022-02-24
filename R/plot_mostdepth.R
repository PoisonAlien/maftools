#'Plot results from mosdepth output
#' @param bed mosdepth output
#' @param col Colors. Default c("#95a5a6", "#7f8c8d")
#' @export

plot_mosdepth = function(bed = NULL, col = c("#95a5a6", "#7f8c8d")){
  
  tum_cov = data.table::fread(input = bed)
  colnames(tum_cov) = c("chr", "start", "end", "doc")
  
  contigs = c(1:22, "X", "Y", paste0("chr", 1:22), "chrX", "chrY")
  tum_cov = tum_cov[chr %in% contigs]
  
  tnxy = tum_cov[chr %in% c('X', 'Y', 'chrX', 'chrY')]
  tn = tum_cov[!chr %in% c('X', 'Y', 'chrX', 'chrY')]
  tn[, chr := gsub(pattern = "chr", replacement = "", x = chr)]
  tn[, chr := as.numeric(as.character(chr))]
  tn = tn[order(chr, start)]
  med_cov = median(tn[,doc], na.rm = TRUE)
  
  
  all_depth = rbind(tn, tnxy)
  colnames(all_depth)[1:3] = c("Chromosome", "Start_Position", "End_Position")
  
  chr.lens.dt = all_depth[,max(End_Position, na.rm = TRUE), .(Chromosome)]
  chr.lens = chr.lens.dt$V1
  names(chr.lens) = chr.lens.dt$Chromosome
  
  #chr.lens = unlist(lapply(split(all_depth, all_depth$Chromosome), function(x) max(x$End_Position, na.rm = TRUE)))
  
  #all_depth = .transformSegments(segmentedData = all_depth, chr.lens = chr.lens)
  
  #log2 median centered
  all_depth[, doc_norm := log2(doc) - log2(med_cov)]
  
  cols = rep(x = col, length(chr.lens))
  
  all_depth_spl = split(all_depth, all_depth$Chromosome)
  all_depth_spl = all_depth_spl[names(chr.lens)]
  
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
    seg.spl.transformed = rbind(seg.spl.transformed, x.seg, fill = TRUE)
  }
  
  all_depth_spl = split(seg.spl.transformed, seg.spl.transformed$Chromosome)
  all_depth_spl = all_depth_spl[names(chr.lens)]
  
  rm(seg.spl.transformed)
  
  cat("Plotting..")
  
  #png(filename = paste0(basename(bed), ".png"), width = 1024, height = 600, bg = "white")
  
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
  })
  abline(v = cumsum(as.numeric(chr.lens)), lty = 2)
  axis(side = 1, at = cumsum(as.numeric(chr.lens)), labels = names(chr.lens))
  axis(side = 2, at = seq(-3, 3, 1), las = 2)
  title(main = "DOC Median centered")
  
  #dev.off()
  
}