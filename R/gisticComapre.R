#' compare two GISTIC objects
#' @details Performs fisher test on 2x2 contingency table generated from two GISTIC objects
#'
#' @param g1 first \code{\link{GISTIC}} object
#' @param g2 second \code{\link{GISTIC}} object
#' @param g1Name optional name for first cohort
#' @param g2Name optional name for second cohort
#' @param minEvent Consider only cytobands with minimum this number of samples altered in at least one of the cohort for analysis. Helpful to ignore single mutated genes. Default 5.
#' @param pseudoCount If TRUE, adds 1 to the contingency table with 0's to avoid `Inf` values in the estimated odds-ratio.
#' @return result list
#' @export
#' @seealso \code{\link{forestPlot}}
#' @seealso \code{\link{lollipopPlot2}}

gisticCompare = function(g1, g2, g1Name = NULL, g2Name = NULL, minEvent = 5, pseudoCount = FALSE){

  g1.gs <- getCytobandSummary(x = g1)
  g2.gs <- getCytobandSummary(x = g2)


  if(is.null(g1Name)){
    g1Name = 'G1'
  }

  if(is.null(g2Name)){
    g2Name = 'G2'
  }

  g1.sampleSize = as.numeric(g1@summary[ID %in% "Samples", summary])
  g2.sampleSize = as.numeric(g2@summary[ID %in% "Samples", summary])

  g1_nums = .getAmpDelCounts(g = g1)
  g2_nums = .getAmpDelCounts(g = g2)

  amptbl = merge(g1_nums[Variant_Classification %in% "Amp", .N , .(Cytoband)], g2_nums[Variant_Classification %in% "Amp", .N , .(Cytoband)], by = "Cytoband", all = TRUE)
  colnames(amptbl) = c("Cytoband", "G1", "G2")
  amptbl[is.na(amptbl)] = 0

  deltbl = merge(g1_nums[Variant_Classification %in% "Del", .N , .(Cytoband)], g2_nums[Variant_Classification %in% "Del", .N , .(Cytoband)], by = "Cytoband", all = TRUE)
  colnames(deltbl) = c("Cytoband", "G1", "G2")
  deltbl[is.na(deltbl)] = 0

  cnvtbl = data.table::rbindlist(l = list(Amp = amptbl, Del = deltbl), idcol = "CNV")
  cnvtbl[,G1_wt := g1.sampleSize - G1]
  cnvtbl[,G2_wt := g2.sampleSize - G2]

  fisherTable = lapply(seq_len(nrow(cnvtbl)), function(i){
    gene = cnvtbl[i, 1]

    ft_mat = matrix(c(cnvtbl[i, G1], cnvtbl[i, G1_wt], cnvtbl[i, G2], cnvtbl[i, G2_wt]),
                    byrow = TRUE, nrow = 2)

    if(length(which(x = ft_mat == 0)) > 0){
      if(pseudoCount){
        ft_mat = ft_mat + 1
      }
    }

    xf = fisher.test(ft_mat, conf.int = TRUE, conf.level = 0.95)

    pval = xf$p.value
    or = xf$estimate
    ci.up = xf$conf.int[2]
    ci.low = xf$conf.int[1]
    tdat = data.table::data.table(pval = pval, or = or, ci.up = ci.up, ci.low = ci.low)
    tdat
  })
  fisherTable = data.table::rbindlist(fisherTable, use.names = TRUE, fill = TRUE)
  fisherTable = cbind(cnvtbl, fisherTable)
  fisherTable = fisherTable[order(pval)]
  fisherTable[,adjPval := p.adjust(p = pval, method = 'fdr')]
  colnames(fisherTable)[3:6] = c(g1Name, g2Name, paste0(g1Name, "_wt"), paste0(g2Name, "_wt"))

  fisherTable
}


.getAmpDelCounts = function(g){
  gs <- getCytobandSummary(x = g)
  amp_bands = gs[Variant_Classification %in% "Amp"]
  del_bands = gs[Variant_Classification %in% "Del"]

  amp_bands_samps = lapply(split(amp_bands, amp_bands$Cytoband), function(cb){

    cb_tsbs = lapply(cb$Unique_Name, function(b){
      b_samps = g@cnMatrix[b,]
      names(b_samps[b_samps != ""])
    })

    cb_tsbs = unique(unlist(cb_tsbs, use.names = FALSE))
    data.table::data.table(Cytoband = unique(cb$Cytoband), Tumor_Sample_Barcode = cb_tsbs)
  })
  amp_bands_samps_tbl = data.table::rbindlist(amp_bands_samps)

  del_bands_samps = lapply(split(del_bands, del_bands$Cytoband), function(cb){

    cb_tsbs = lapply(cb$Unique_Name, function(b){
      b_samps = g@cnMatrix[b,]
      names(b_samps[b_samps != ""])
    })

    cb_tsbs = unique(unlist(cb_tsbs, use.names = FALSE))
    data.table::data.table(Cytoband = unique(cb$Cytoband), Tumor_Sample_Barcode = cb_tsbs)
  })
  del_bands_samps_tbl = data.table::rbindlist(del_bands_samps)

  data.table::rbindlist(list(Amp = amp_bands_samps_tbl, Del= del_bands_samps_tbl), idcol = "Variant_Classification")
}
