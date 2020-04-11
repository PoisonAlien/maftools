#' Compare genes to known TCGA drivers and their biological pathways
#' @description A small function which uses known cancer driver genes and their associatd pathways from TCGA cohorts. See reference for details
#' @param m an \code{\link{MAF}} object
#' @param top Top number of genes to use. Default 50
#' @references Bailey MH, Tokheim C, Porta-Pardo E, et al. Comprehensive Characterization of Cancer Driver Genes and Mutations . Cell. 2018;173(2):371–385.e18. doi:10.1016/j.cell.2018.02.060
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' tcgaDriverBP(m = laml)
#' @export
tcgaDriverBP = function(m, top = 50){
  pathways = system.file("extdata", "BP_SMGs.txt.gz", package = "maftools")
  pathways = data.table::fread(file = pathways, skip = "Gene")
  pathways = split(pathways, as.factor(as.character(pathways$Gene)))
  pathways = lapply(pathways, function(x){
    tcga_driver = paste(x$Cancer_type, collapse = ",")
    x = x[!duplicated(Gene)]
    x[,tcga_driver := tcga_driver]
    x
  })
  pathways = data.table::rbindlist(l = pathways)
  pathways[, Cancer_type := NULL]
  pathways$Pathway = gsub(pattern = " ", replacement = "_", x = pathways$Pathway)

  temp_dat = getGeneSummary(x = m)[1:top,.(Hugo_Symbol, AlteredSamples)]
  temp_dat = merge(temp_dat, pathways, all.x = TRUE, by.x = "Hugo_Symbol", by.y = "Gene")
  temp_dat[is.na(temp_dat)] = "Unknown"
  temp_dat = temp_dat[order(-AlteredSamples, -Pathway, Hugo_Symbol)]
  temp_dat = rbind(temp_dat[!Pathway %in% "Unknown"], temp_dat[Pathway %in% "Unknown"])

  nsamps = as.numeric(m@summary[ID == 'Samples', summary])

  temp_dat[, pctAltered := paste0(round(AlteredSamples/nsamps * 100, 2), "%")]
  temp_dat = temp_dat[,.(Hugo_Symbol, AlteredSamples, pctAltered, Pathway, tcga_driver)]
  temp_dat
}
