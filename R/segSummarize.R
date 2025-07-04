#' Summarize CBS segmentation results
#' @details
#' A handy function to summarize CBS segmentation results. Takes segmentation results generated by DNAcopy package \code{\link{segment}} and summarizes the CN for each cytoband, chromosomal arms, and whole chromosomes.
#'
#' @param seg segmentation results generated from \code{DNAcopy} package \code{\link{segment}}. Input should be a multi-sample segmentation file or a data.frame. First six columns should correspond to sample name, chromosome,  start,  end, Num_Probes, Segment_Mean in log2 scale. (default output format from DNAcopy)
#' @param build genome build. Default hg19. Can be hg19, hg38. If other than these, use `cytoband` argument
#' @param cytoband cytoband data from UCSC genome browser. Only needed if `build` is other than `hg19` or `hg38`
#' @param thr threshold to call amplification and deletion. Any cytobands or chromosomal arms with median logR above or below this will be called. Default 0.3
#' @param chr_coverage minimum percentage of chromosome that must exceed threshold to call whole chromosome gain/loss. Default 0.75 (75%)
#' @param filter_arm_events logical. If TRUE (default), removes arm-level events where both p and q arms are affected in the same direction, reclassifying them as chromosome-level events. This prevents double-counting and follows standard cytogenetic classification. Set to FALSE to retain the original behavior where each arm is counted independently.
#' @param verbose Default TRUE
#' @param maf optional MAF
#' @param plot Default TRUE. Draw heatmap
#' @param genes Add mutation status of these genes as an annotation to the heatmap
#' @param topanno annotation for each sample. This is passed as an input to `annotation_col` of `pheatmap`
#' @param topannocols annotation cols for `topanno`. This is passed as an input to `annotation_colors` of `pheatmap`
#' @return List of median CN values for each cytoband, chromosomal arm, and whole chromosome along with the plotting matrix
#' @export
#' @examples
#' laml.seg <- system.file("extdata", "LAML_CBS_segments.tsv.gz", package = "maftools")
#' segSummarize(seg = laml.seg)
#'
#' #Heighlight some genes as annotation
#' laml.maf = system.file("extdata", "tcga_laml.maf.gz", package = "maftools") #MAF file
#' laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') #clinical data
#' laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
#'
#' segSummarize(seg = laml.seg, maf = laml, genes = c("FLT3", "DNMT3A"))
#'
#' #Use stricter criteria for whole chromosome events
#' segSummarize(seg = laml.seg, chr_coverage = 0.9)
#'
#' #Disable arm event filtering to use original behavior
#' segSummarize(seg = laml.seg, filter_arm_events = FALSE)

segSummarize = function(seg = NULL, build = "hg19", cytoband = NULL, thr = 0.3, chr_coverage = 0.75, filter_arm_events = TRUE, verbose = TRUE, maf = NULL, genes = NULL, topanno = NULL, topannocols = NA, plot = TRUE){

  if(is.null(cytoband)){
    build = match.arg(arg = build, choices = c("hg19", "hg38"))
    cytoband <- system.file("extdata", paste0(build, "_cytobands.tsv.gz"), package = "maftools")
    if(!file.exists(cytoband)){
      stop("Cytoband file does not exist!")
    }
  }

  if (is(object = seg, class2 = "data.frame")) {
    seg = data.table::as.data.table(seg)
  }else {
    seg = data.table::fread(input = seg)
  }
  colnames(seg)[1:6] = c("Sample", "Chromosome",  "Start",  "End", "Num_Probes", "Segment_Mean")
  seg$Chromosome = gsub(pattern = "^chr", replacement = "", x = seg$Chromosome)
  data.table::setkey(x = seg, Chromosome, Start, End)

  if(nrow(seg[,.N,Sample]) < 3){
    stop("Not enough samples to proceed (N: ", nrow(seg[,.N,Sample]), ")")
  }

  cy = data.table::fread(input = cytoband)
  colnames(cy) = c("Chromosome", "Start", "End", "name", "stain")
  cy$arm = substr(cy$name, 1, 1)
  cy$Chromosome = gsub(pattern = "^chr", replacement = "", x = cy$Chromosome)
  data.table::setkey(x = cy, Chromosome, Start, End)

  seg_events = lapply(split(seg, seg$Sample), function(x){
    .seg2arm(se = x, cy = cy, thr = thr, chr_coverage = chr_coverage, filter_arm_events = filter_arm_events)
  })

  # Collect reclassification information for verbose output
  if (filter_arm_events & verbose) {
    reclassified_events = lapply(names(seg_events), function(sample_name) {
      reclassified = seg_events[[sample_name]]$reclassified_info
      if (!is.null(reclassified) && nrow(reclassified) > 0) {
        reclassified[, sample := sample_name]
        return(reclassified)
      }
      return(NULL)
    })
    reclassified_events = reclassified_events[!sapply(reclassified_events, is.null)]

    if (length(reclassified_events) > 0) {
      all_reclassified = data.table::rbindlist(reclassified_events)
      reclassified_summary = all_reclassified[, .N, .(chromosome, classification)]
      message("Reclassified arm events (both p+q arms affected): ",
              paste(paste0("chr", reclassified_summary$chromosome, "(", reclassified_summary$classification, "):", reclassified_summary$N), collapse = ", "))
    }
  }

  cy[, id := paste(Chromosome, name, sep = "_")]
  cy[, arm := paste(Chromosome, substr(name, 1, 1), sep = "_")]

  #Summarize CN per chromosomal cytoband
  cytoband_events = lapply(seg_events, function(x) x$seg) |> data.table::rbindlist(idcol = "Tumor_Sample_Barcode")
  cytoband_events = data.table::dcast(data = cytoband_events, formula = chromosome+arm ~ Tumor_Sample_Barcode, value.var = "cn")

  cytoband_events[, id := paste(chromosome, arm, sep = "_")]
  cytoband_events$chromosome = NULL
  cytoband_events$arm = NULL
  data.table::setDF(x = cytoband_events, cytoband_events$id)
  cytoband_events$id = NULL
  cytoband_events = cytoband_events[cy[id %in% rownames(cytoband_events), id],,]
  cytoband_events[cytoband_events > 4] = 4
  cytoband_events = cytoband_events[complete.cases(cytoband_events),]

  # Note: Chromosome-level events are kept separate and returned in chromosome_events
  # The heatmap matrix focuses on cytoband-level resolution only

  ##row anno
  cn_band_anno = as.data.frame(cy[id %in% rownames(cytoband_events)])
  cn_band_anno = data.frame(row.names = cn_band_anno$id, chr = cn_band_anno$Chromosome, arm = substr(x = cn_band_anno$name, 1, 1))

  # No chromosome-level annotations needed since we focus on cytoband resolution

  #order by chromosome for human build
  chr_ord = intersect(c(1:22, "X", "Y"), names(split(cn_band_anno, cn_band_anno$chr)))
  cn_band_anno = cn_band_anno[unlist(lapply(split(cn_band_anno, cn_band_anno$chr)[chr_ord], rownames), use.names = FALSE), , drop = FALSE]
  cytoband_events = cytoband_events[rownames(cn_band_anno),, drop = FALSE]

  ##Heatmap colors
  cols_chromosome = setNames(object = rep(c("#95a5a6", "#7f8c8d"), length(unique(cn_band_anno$chr))), nm = unique(cn_band_anno$chr))
  cols_chromosome = cols_chromosome[unique(cn_band_anno$chr)]
  cols_chromosome_arms = setNames(object = c("#2c3e50", "#ecf0f1", "#bdc3c7"), nm = c("whole", "p", "q"))

  gene2barcode = gene2barcode_cols = NA
  show_annotation_legend = FALSE

  if(!is.null(maf)){
    if(!is.null(genes)){
      gene2barcode = genesToBarcodes(maf = maf, genes = genes, justNames = TRUE)
      print(gene2barcode)
      gene2barcode = lapply(gene2barcode, function(x) data.table::data.table(Tumor_Sample_Barcode = x, status = 'yes')) |> data.table::rbindlist(idcol = "Hugo_Symbol") |> data.table::dcast(formula = Tumor_Sample_Barcode ~ Hugo_Symbol, value.var = "status")
      data.table::setDF(x = gene2barcode, rownames = gene2barcode$Tumor_Sample_Barcode)
      gene2barcode = gene2barcode[, 2:ncol(gene2barcode), drop = FALSE]
      gene2barcode_cols = lapply(X = colnames(gene2barcode), function(x) setNames(object = "#34495e", nm = "yes"))
      names(gene2barcode_cols) = colnames(gene2barcode)
    }
  }

  if(!is.null(topanno)){
    if(!is.null(nrow(gene2barcode))){
      gene2barcode = merge(gene2barcode, topanno, by = "row.names", all = TRUE)
      rownames(gene2barcode) = gene2barcode$`Row.names`
      gene2barcode$`Row.names` = NULL
    }else{
      gene2barcode = topanno
    }
    show_annotation_legend = TRUE
  }

  if(plot){
    pheatmap::pheatmap(
      cytoband_events,
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      show_colnames = FALSE,
      annotation_row = cn_band_anno,
      color = grDevices::colorRampPalette(c(
        "#08519C",  # CN=0: Dark blue (deep loss)
        "#3182BD",  # CN=1: Medium blue (loss)
        "#FFFFFF",  # CN=2: White (normal diploid)
        "#DE2D26",  # CN=3: Medium red (gain)
        "#A50F15"   # CN=4+: Dark red (high amplification)
      ))(100),
      annotation_colors = c(list(
        chr = cols_chromosome,
        arm = cols_chromosome_arms
      ), gene2barcode_cols, topannocols),
      show_rownames = FALSE,
      border_color = 'black', annotation_legend = show_annotation_legend, legend_breaks = seq(0, 4, 1), legend_labels = c("0", "1", "2", "3", ">4"), annotation_col = gene2barcode, na_col = "gray", breaks = seq(0, 4, length.out = 101)
    )
  }

  #Arm events
  arm_events = lapply(seg_events, function(x) x$arm) |> data.table::rbindlist(fill = TRUE, use.names = TRUE, idcol = "Tumor_Sample_Barcode")
  arm_events = arm_events[!is.na(Variant_Classification)][order(chromosome)]
  arm_events_n = arm_events[,.N,.(arm,Variant_Classification)][!is.na(Variant_Classification)][order(-N)]

  if(verbose){
    if(filter_arm_events){
      message("Recurrent chromosomal arm aberrations (filtered: arms with concordant p/q events reclassified as chromosome-level)")
    } else {
      message("Recurrent chromosomal arm aberrations (unfiltered: each arm counted independently)")
    }
    print(arm_events_n)
  }

  #Focal band events
  band_events = lapply(seg_events, function(x) x$seg) |> data.table::rbindlist(fill = TRUE, use.names = TRUE, idcol = "Tumor_Sample_Barcode")
  band_events = band_events[!is.na(Variant_Classification)][order(chromosome)]
  colnames(band_events)[3] = "cytoband"
  band_events_n = band_events[,.N,.(chromosome, cytoband,Variant_Classification)][!is.na(Variant_Classification)][order(-N)]

  #Chromosome-level events
  chromosome_events = lapply(seg_events, function(x) x$chromosome) |> data.table::rbindlist(fill = TRUE, use.names = TRUE, idcol = "Tumor_Sample_Barcode")
  chromosome_events = chromosome_events[!is.na(Variant_Classification)][order(chromosome)]
  chromosome_events_n = chromosome_events[,.N,.(chromosome,Variant_Classification)][!is.na(Variant_Classification)][order(-N)]

  if(verbose){
    message("Recurrent whole chromosome aberrations")
    print(chromosome_events_n)
  }

  attr(cytoband_events, "meta") = list(top_anno = gene2barcode, top_anno_cols = gene2barcode_cols)

  list(heatmap_matrix = cytoband_events, arm_events = arm_events, band_events = band_events, chromosome_events = chromosome_events)
}

.seg2arm = function (se, cy = NULL, thr = 0.3, chr_coverage = 0.75, filter_arm_events = TRUE) {

  data.table::setkey(x = se, Chromosome, Start, End)
  se_olaps = data.table::foverlaps(x = cy, y = se)
  cytoband_mean = se_olaps[, median(Segment_Mean, na.rm = TRUE),
                           .(Chromosome, name)]
  colnames(cytoband_mean) = c("chromosome", "arm", "logR")
  cytoband_mean[, `:=`(cn, 2 * (2^logR))]
  cytoband_mean$Variant_Classification = ifelse(cytoband_mean$logR >
                                                  thr, "Amp", no = ifelse(cytoband_mean$logR < -thr, yes = "Del",
                                                                          no = NA))
  cytoband_mean = cytoband_mean[!chromosome %in% c("chrX",
                                                   "chrY")]
  arm_mean = se_olaps[, median(Segment_Mean, na.rm = TRUE),
                      .(Chromosome, arm)]
  colnames(arm_mean) = c("chromosome", "arm", "logR")
  arm_mean[, `:=`(cn, 2 * (2^logR))]
  arm_mean$Variant_Classification = ifelse(arm_mean$logR >
                                             thr, "Gain", no = ifelse(arm_mean$logR < -thr, yes = "Loss",
                                                                     no = NA))
  arm_mean$arm = paste(arm_mean$chromosome, arm_mean$arm, sep = "_")
  arm_mean = arm_mean[!chromosome %in% c("chrX", "chrY")]

  # Calculate chromosome-level events using percentage coverage approach
  chr_events = .seg2chromosome(se = se, cy = cy, thr = thr, chr_coverage = chr_coverage)

  # Filter arm events if both arms are affected in same direction
  reclassified_info = data.table::data.table()
  if (filter_arm_events) {
    # Identify chromosomes where both arms are affected in same direction
    arm_concordance = .detect_arm_concordance(arm_mean, thr = thr)

    # Filter out concordant arm events (they become chromosome events)
    arm_mean_filtered = arm_mean[!arm %in% arm_concordance$concordant_arms]

    # Add reclassified arm events to chromosome events
    reclassified_chr_events = arm_concordance$reclassified_chr_events
    if (nrow(reclassified_chr_events) > 0) {
      # Ensure column structure matches chr_events
      if ("Chromosome" %in% names(chr_events)) {
        reclassified_chr_events[, Chromosome := chromosome]
      }
      chr_events = rbind(chr_events, reclassified_chr_events, fill = TRUE)
    }

    arm_mean = arm_mean_filtered
    reclassified_info = arm_concordance$reclassified_info
  }

  list(arm = arm_mean, seg = cytoband_mean, chromosome = chr_events, reclassified_info = reclassified_info)
}

.seg2chromosome = function(se, cy = NULL, thr = 0.3, chr_coverage = 0.75) {

  data.table::setkey(x = se, Chromosome, Start, End)

  # Calculate chromosome lengths for coverage calculation
  chr_lengths = se[, .(chr_length = max(End) - min(Start) + 1), by = Chromosome]

  # For each chromosome, calculate percentage of segments above/below threshold
  chr_events = se[, {
    total_length = max(End) - min(Start) + 1

    # Calculate length-weighted segments above/below threshold
    gain_length = sum(ifelse(Segment_Mean > thr, End - Start + 1, 0))
    loss_length = sum(ifelse(Segment_Mean < -thr, End - Start + 1, 0))

    # Calculate coverage percentages
    gain_coverage = gain_length / total_length
    loss_coverage = loss_length / total_length

    # Calculate median logR for final classification
    median_logR = median(Segment_Mean, na.rm = TRUE)

    # Determine classification based on coverage and median
    classification = NA_character_
    if (gain_coverage >= chr_coverage & median_logR > thr) {
      classification = "Gain"
    } else if (loss_coverage >= chr_coverage & median_logR < -thr) {
      classification = "Loss"
    }

    # Calculate copy number
    cn = 2 * (2^median_logR)

    list(
      chromosome = Chromosome[1],
      logR = median_logR,
      cn = cn,
      gain_coverage = gain_coverage,
      loss_coverage = loss_coverage,
      Variant_Classification = classification
    )
  }, by = Chromosome]

  # Check arm concordance - both p and q arms must be in same direction
  if (!is.null(cy)) {
    data.table::setkey(x = cy, Chromosome, Start, End)
    cy_arms = cy[, .(arm = unique(arm)), by = Chromosome]

    # Get arm-level data for concordance check
    se_olaps = data.table::foverlaps(x = cy, y = se)
    arm_summary = se_olaps[, .(arm_logR = median(Segment_Mean, na.rm = TRUE)), .(Chromosome, arm)]

    # Check if both arms have same direction for each chromosome
    arm_concordance = arm_summary[, {
      p_arm_logR = arm_logR[arm == "p"]
      q_arm_logR = arm_logR[arm == "q"]

      # Both arms must be in same direction and exceed threshold
      concordant = FALSE
      if (length(p_arm_logR) > 0 & length(q_arm_logR) > 0) {
        concordant = (p_arm_logR > thr & q_arm_logR > thr) | (p_arm_logR < -thr & q_arm_logR < -thr)
      }

      list(arm_concordant = concordant)
    }, by = Chromosome]

    # Merge concordance info and filter
    chr_events = merge(chr_events, arm_concordance, by = "Chromosome", all.x = TRUE)
    chr_events[is.na(arm_concordant) | !arm_concordant, Variant_Classification := NA]
    chr_events$arm_concordant = NULL
  }

  # Exclude sex chromosomes and filter out neutral events
  chr_events = chr_events[!chromosome %in% c("chrX", "chrY") & !is.na(Variant_Classification)]

  return(chr_events)
}

.detect_arm_concordance = function(arm_mean, thr = 0.3) {
  # arm_mean has columns: chromosome, arm, logR, cn, Variant_Classification
  # arm column format: "1_p", "1_q", etc.

  # Extract chromosome and arm info
  arm_dt = data.table::as.data.table(arm_mean)
  arm_dt[, c("chr_num", "arm_type") := data.table::tstrsplit(arm, "_", fixed = TRUE)]

  # Only consider arms with significant events (not NA)
  arm_events = arm_dt[!is.na(Variant_Classification)]

  if (nrow(arm_events) == 0) {
    return(list(concordant_arms = character(0), reclassified_chr_events = data.table::data.table()))
  }

  # Check for concordance within each chromosome
  concordant_chrs = arm_events[, {
    p_events = .SD[arm_type == "p"]
    q_events = .SD[arm_type == "q"]

    # Check if both arms exist and have same direction
    if (nrow(p_events) > 0 & nrow(q_events) > 0) {
      p_classification = p_events$Variant_Classification[1]
      q_classification = q_events$Variant_Classification[1]

      # Both arms must be in same direction (both Gain or both Loss)
      if (p_classification == q_classification) {
        list(
          concordant = TRUE,
          classification = p_classification,
          median_logR = median(c(p_events$logR[1], q_events$logR[1])),
          p_arm = p_events$arm[1],
          q_arm = q_events$arm[1]
        )
      } else {
        list(concordant = FALSE, classification = NA_character_, median_logR = NA_real_, p_arm = NA_character_, q_arm = NA_character_)
      }
    } else {
      list(concordant = FALSE, classification = NA_character_, median_logR = NA_real_, p_arm = NA_character_, q_arm = NA_character_)
    }
  }, by = chr_num]

  # Get concordant arms to remove
  concordant_arms = character(0)
  reclassified_chr_events = data.table::data.table()
  concordant_subset = concordant_chrs[concordant == TRUE]

  if (nrow(concordant_subset) > 0) {
    # Arms to remove from arm-level analysis
    concordant_arms = c(concordant_subset$p_arm, concordant_subset$q_arm)

    # Create chromosome-level events from reclassified arms
    reclassified_chr_events = concordant_subset[, .(
      chromosome = chr_num,
      logR = median_logR,
      cn = 2 * (2^median_logR),
      gain_coverage = ifelse(classification == "Gain", 1.0, 0.0),
      loss_coverage = ifelse(classification == "Loss", 1.0, 0.0),
      Variant_Classification = classification
    )]
  }

  return(list(
    concordant_arms = concordant_arms,
    reclassified_chr_events = reclassified_chr_events,
    reclassified_info = if(nrow(concordant_subset) > 0) concordant_subset[, .(chromosome = chr_num, classification)] else data.table::data.table()
  ))
}
