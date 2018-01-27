#' compare two cohorts (MAF).
#' @details Performs fisher test on 2x2 contigency table generated from two cohorts to find differentially mutated genes.
#'
#' @param m1 first \code{\link{MAF}} object
#' @param m2 second \code{\link{MAF}} object
#' @param m1Name optional name for first cohort
#' @param m2Name optional name for second cohort
#' @param minMut Consider only genes with minimum this number of samples mutated in atleast one of the cohort for analysis. Helful to ignore single mutated genes. Default 5.
#' @param useCNV whether to include copy number events to compare MAFs. Only applicable when MAF is read along with copy number data. Default TRUE if available.
#' @return result list
#' @examples
#' primary.apl <- system.file("extdata", "APL_primary.maf.gz", package = "maftools")
#' relapse.apl <- system.file("extdata", "APL_relapse.maf.gz", package = "maftools")
#' primary.apl <- read.maf(maf = primary.apl)
#' relapse.apl <- read.maf(maf = relapse.apl)
#' pt.vs.rt <- mafCompare(m1 = primary.apl, m2 = relapse.apl, m1Name = 'Primary',
#' m2Name = 'Relapse', minMut = 5)
#' @export
#' @seealso \code{\link{forestPlot}}

mafCompare = function(m1, m2, m1Name = NULL, m2Name = NULL, minMut = 5, useCNV = TRUE){

  m1.gs <- getGeneSummary(x = m1)
  m2.gs <- getGeneSummary(x = m2)


   if(is.null(m1Name)){
     m1Name = 'M1'
   }

   if(is.null(m2Name)){
     m2Name = 'M2'
   }

  if(useCNV){
    m1.genes = as.character(m1.gs[AlteredSamples >= minMut,Hugo_Symbol])
    m2.genes = as.character(m2.gs[AlteredSamples >= minMut,Hugo_Symbol])
    uniqueGenes = unique(c(m1.genes, m2.genes))
  }else{
    m1.genes = as.character(m1.gs[MutatedSamples >= minMut, Hugo_Symbol])
    m2.genes = as.character(m2.gs[MutatedSamples >= minMut, Hugo_Symbol])
    uniqueGenes = unique(c(m1.genes, m2.genes))
  }

 #com.genes = intersect(m1.gs[,Hugo_Symbol], m2.gs[,Hugo_Symbol])

 m1.sampleSize = as.numeric(m1@summary[3, summary])
 m2.sampleSize = as.numeric(m2@summary[3, summary])

 m1.gs.comGenes = m1.gs[Hugo_Symbol %in% uniqueGenes]
 m2.gs.comGenes = m2.gs[Hugo_Symbol %in% uniqueGenes]

 sampleSummary = data.table::data.table(Cohort = c(m1Name, m2Name), SampleSize = c(m1.sampleSize, m2.sampleSize))

 if(useCNV){
   m.gs.meged = merge(m1.gs.comGenes[,.(Hugo_Symbol, AlteredSamples)], m2.gs.comGenes[,.(Hugo_Symbol, AlteredSamples)],
                      by = 'Hugo_Symbol', all = TRUE)
 }else{
   m.gs.meged = merge(m1.gs.comGenes[,.(Hugo_Symbol, MutatedSamples)], m2.gs.comGenes[,.(Hugo_Symbol, MutatedSamples)],
                      by = 'Hugo_Symbol', all = TRUE)
 }

 #Set missing genes to zero
 m.gs.meged[is.na(m.gs.meged)] = 0

 m.gs.meged = as.data.frame(m.gs.meged)

 fisherTable = c()

 for(i in 1:nrow(m.gs.meged)){
   gene = m.gs.meged[i, 1]
   m1Mut = m.gs.meged[i,2]
   m2Mut = m.gs.meged[i,3]
   #print(i)
   xf = fisher.test(matrix(c(m1Mut, m1.sampleSize-m1Mut, m2Mut, m2.sampleSize-m2Mut),
                           byrow = TRUE, nrow = 2), conf.int = TRUE, conf.level = 0.95)

   pval = xf$p.value
   or = xf$estimate
   ci.up = xf$conf.int[1]
   ci.low = xf$conf.int[2]
   tdat = data.table::data.table(Hugo_Symbol = gene, m1Mut , m2Mut, pval = pval, or = or, ci.up = ci.up, ci.low = ci.low)
   fisherTable = rbind(fisherTable, tdat)
 }

  fisherTable = fisherTable[order(pval)]
  fisherTable[,adjPval := p.adjust(p = pval, method = 'fdr')]

  colnames(fisherTable)[2:3] = c(m1Name, m2Name)

 return(list(results = fisherTable, SampleSummary = sampleSummary))

}
