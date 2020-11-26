#' Performs mutational enrichment analysis for a given clinical feature.
#' @description Performs pairwise and groupwise fisher exact tests to find differentially enriched genes for every factor within a clinical feature.
#'
#' @param maf \code{\link{MAF}} object
#' @param clinicalFeature columns names from `clinical.data` slot of \code{MAF} to be analysed for.
#' @param minMut Consider only genes with minimum this number of samples mutated. Default 5.
#' @param useCNV whether to include copy number events. Only applicable when MAF is read along with copy number data. Default TRUE if available.
#' @param annotationDat If MAF file was read without clinical data, provide a custom \code{data.frame} or a tsv file with a column containing Tumor_Sample_Barcodes along with clinical features. Default NULL.
#' @return result list containing p-values
#' @details Performs fishers test on 2x2 contingency table for WT/Mutants in group of interest vs rest of the sample. Odds Ratio indicate the odds of observing mutant in the group of interest compared to wild-type
#' @examples
#' \dontrun{
#' laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
#' laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')
#' laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
#' clinicalEnrichment(laml, 'FAB_classification')
#' }
#' @seealso \code{\link{plotEnrichmentResults}}
#' @export

clinicalEnrichment = function(maf, clinicalFeature = NULL, annotationDat = NULL, minMut = 5, useCNV = TRUE){

  if(is.null(clinicalFeature)){
    stop("Missing clinicalFeature. Use getClinicalData() to see available features.")
  }

  if(is.null(annotationDat)){
    cd = getClinicalData(x = maf)[,c("Tumor_Sample_Barcode", clinicalFeature), with = FALSE]
  }else{
    if(is.data.frame(annotationDat)){
      cd = data.table::as.data.table(annotationDat)
      cd = cd[,c("Tumor_Sample_Barcode", clinicalFeature), with = FALSE]
    }else if(file.exists(annotationDat)){
      cd = data.table::fread(input = annotationDat)
      cd = cd[,c("Tumor_Sample_Barcode", clinicalFeature), with = FALSE]
    }
  }


  colnames(cd)[2] = 'cf'
  cd$cf = as.character(cd$cf)
  cf.tbl = table(cd$cf)
  message(paste0("Sample size per factor in ", clinicalFeature, ":"))
  print(cf.tbl)

  if(length(cf.tbl) == 1){
    stop("Single factor. Nothing to compare..")
  }

  #Source code from reporttools (https://github.com/cran/reporttools/blob/master/R/pairwise.fisher.test.r)
  pairwise.fisher.test <- function(x, g, p.adjust.method, ...){
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    g <- factor(g)

    compare.levels <- function(i, j) {
      ind <- (as.integer(g) %in% c(i, j)) & (is.na(x) == FALSE) & (is.na(g) == FALSE)
      xi <- factor(x[ind], exclude = NULL)
      xj <- factor(g[ind], exclude = NULL)
      tab <- table(xi, xj)
      nonzeromarginal <- (min(apply(tab, 1, sum)) * min(apply(tab, 2, sum)) > 0)
      size <- ((nrow(tab) > 1) * (ncol(tab) > 1) > 0)
      if ((nonzeromarginal == TRUE) & (size == TRUE)){fisher.test(xi, xj, ...)$p.value} else {NA}
    }

    PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
    ans <- list(data.name = DNAME, p.value = PVAL, p.adjust.method = p.adjust.method)
    class(ans) <- "pairwise.htest"
    return(ans)
  }

  if(useCNV){
    genes = as.character(getGeneSummary(x = maf)[AlteredSamples > minMut, Hugo_Symbol])
  }else{
    genes = as.character(getGeneSummary(x = maf)[MutatedSamples > minMut, Hugo_Symbol])
  }

  plist = lapply(genes, function(x){
          g = unique(genesToBarcodes(maf = maf, genes = x, justNames = TRUE)[[1]])
          cd$Genotype = ifelse(test = cd$Tumor_Sample_Barcode %in% g, yes = "Mutant", no = "WT")

          #Perform groupwise comparison for given gene
          ft = lapply(X = names(cf.tbl), FUN = function(y){
            cd$Group = ifelse(test = cd$cf %in% y, yes = y, no = "Other")
            cd$Genotype = factor(x = cd$Genotype, levels = c("WT", "Mutant"))
            cd$Group = factor(x = cd$Group, levels = c(y, "Other"))
            cd.tbl = with(cd, table(Genotype, Group))
            cd.tbl = cd.tbl[c("Mutant", "WT") ,c(y, "Other")]
            ft = fisher.test(cd.tbl)
            ft.tbl = data.table::data.table(Group1 = y, Group2 = "Rest",
                                            n_mutated_group1 = paste0(nrow(cd[Group %in% y][Genotype %in% 'Mutant']), " of ", nrow(cd[Group %in% y])),
                                            n_mutated_group2 = paste0(nrow(cd[!Group %in% y][Genotype %in% 'Mutant']), " of ", nrow(cd[!Group %in% y])),
                                            p_value = ft$p.value, OR = ft$estimate, OR_low = ft$conf.int[1], OR_high = ft$conf.int[2],
                                            Hugo_Symbol = x, Analysis = "Group")
            ft.tbl
          })
          ft = data.table::rbindlist(ft)

          #Perform pairwise fisher test for every gene
          prop.tbl = pairwise.fisher.test(x = cd$Genotype, g = cd$cf, p.adjust.method = "fdr")
          ptbl = data.table::as.data.table(as.data.frame(prop.tbl$p.value), keep.rownames = TRUE)
          ptbl = data.table::melt(ptbl, id.vars = "rn")
          colnames(ptbl) = c("Var1", "Var2", "value")
          ptbl[,Hugo_Symbol := x][,Analysis := "Pairwise"]
          ptbl = ptbl[,.(Hugo_Symbol, Var1, Var2, value, Analysis)]
          colnames(ptbl) = c("Hugo_Symbol", "Feature_1", "Feature_2", "fdr", "Analysis")
          ptbl = ptbl[!is.na(fdr)]

          f1.mutants = cd[,.N,.(cf, Genotype)][Genotype %in% 'Mutant', .(cf, N)]
          if(length(names(cf.tbl)[!names(cf.tbl) %in% f1.mutants[,cf]]) > 0){
            #Add zero counts for missing factors
            f1.mutants = rbind(f1.mutants,
                               data.table::data.table(cf = names(cf.tbl)[!names(cf.tbl) %in% f1.mutants[,cf]], N = 0))
          }

          f1.mutants = merge(f1.mutants, cd[,.N,.(cf)], by = 'cf')
          f1.mutants[,n_mutated_Feature := paste0(N.x, " of ", N.y)]

          ptbl = merge(ptbl, f1.mutants[,.(cf, n_mutated_Feature)], by.x = 'Feature_1', by.y = 'cf', all.x = TRUE)
          ptbl = merge(ptbl, f1.mutants[,.(cf, n_mutated_Feature)], by.x = 'Feature_2', by.y = 'cf', all.x = TRUE)
          colnames(ptbl)[6:7] = c('n_mutated_Feature1', 'n_mutated_Feature2')
          ptbl = ptbl[,.(Hugo_Symbol, Feature_1, Feature_2, n_mutated_Feature1, n_mutated_Feature2, fdr, Analysis)]

          ptbl = rbind(ptbl, ft, fill = TRUE)
          ptbl
        })

  plist = data.table::rbindlist(l = plist, fill = TRUE)

  pw.pvals = plist[Analysis %in% "Pairwise",.(Hugo_Symbol, Feature_1, Feature_2, n_mutated_Feature1, n_mutated_Feature2, fdr)][order(fdr)]
  gw.pvals = plist[Analysis %in% "Group",.(Hugo_Symbol, Group1, Group2, n_mutated_group1, n_mutated_group2, p_value, OR, OR_low, OR_high)][order(p_value)]
  gw.pvals[,fdr := p.adjust(p_value, method = "fdr")]

  return(list(pairwise_comparision = pw.pvals, groupwise_comparision = gw.pvals, cf_sizes = cd[,.N,cf], clinicalFeature = clinicalFeature))
}
