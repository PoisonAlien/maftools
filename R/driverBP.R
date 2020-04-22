#' Compare genes to known TCGA drivers and their biological pathways
#' @description A small function which uses known cancer driver genes and their associatd pathways from TCGA cohorts. See reference for details
#' @param m an \code{\link{MAF}} object
#' @param genes genes to compare. Default `NULL`.
#' @param top Top number of genes to use. Mutually exclusive with `genes` argument. Default 20
#' @fontSize Default 0.7
#' @references Bailey MH, Tokheim C, Porta-Pardo E, et al. Comprehensive Characterization of Cancer Driver Genes and Mutations . Cell. 2018;173(2):371–385.e18. doi:10.1016/j.cell.2018.02.060
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' tcgaDriverBP(m = laml)
#' @export
tcgaDriverBP = function(m, genes = NULL, top = 20, fontSize = 0.7){
  pathways = system.file("extdata", "BP_SMGs.txt.gz", package = "maftools")
  pathways = data.table::fread(file = pathways, skip = "Gene")


  if(!is.null(genes)){
    temp_dat = getGeneSummary(x = m)[Hugo_Symbol %in% genes,.(Hugo_Symbol, AlteredSamples)]
  }else{
    temp_dat = getGeneSummary(x = m)[1:top,.(Hugo_Symbol, AlteredSamples)]
  }



  temp_dat = merge(temp_dat, pathways, all.x = TRUE, by.x = "Hugo_Symbol", by.y = "Gene")
  temp_dat[is.na(temp_dat)] = "None"
  temp_cast = data.table::dcast(data = temp_dat, Hugo_Symbol  ~ Cancer_type, value.var = 'Pathway')
  data.table::setDF(x = temp_cast, rownames = temp_cast$Hugo_Symbol)
  temp_cast = temp_cast[,-1]
  temp_cast$None = NULL
  temp_cast[!is.na(temp_cast)] = 1
  temp_cast[is.na(temp_cast)] = 0
  nm = temp_cast[names(sort(apply(temp_cast, 1, function(x) sum(as.numeric(x))), decreasing = TRUE)),,]
  nm[nm == 0] = NA
  nm = t(apply(nm, 2, rev))

  par(mar = c(4, 4, 2, 2), xpd = TRUE)
  image(x = 1:nrow(nm), y = 1:ncol(nm), z = apply(nm, 2, as.numeric), axes = FALSE, xaxt="n", yaxt="n", xlab="", ylab="", col = "brown", ) #col = "#FC8D62"
  abline(h = (1:ncol(nm)) + 0.5, col = "white")
  abline(v = (1:nrow(nm)) + 0.5, col = "white")
  points(which(is.na(nm), arr.ind = TRUE), pch=".", col = "gray70")
  mtext(text = colnames(nm), side = 2, at = 1:ncol(nm),
        font = 3, line = 0.4, cex = fontSize, las = 2)

  text(x =1:nrow(nm), y = par("usr")[3] - 0.2,
       labels = rownames(nm), srt = 90, font = 3, cex = fontSize, adj = 1)
  title(main = "TCGA cohorts", adj = 1)
  #legend(x = "bottomright", legend = "Driver gene", col = "brown", pch = 16, xpd = TRUE)

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

  temp_dat = temp_dat[!duplicated(Hugo_Symbol)]
  temp_dat[,Cancer_type := NULL]
  temp_dat[,Pathway := NULL]
  temp_dat = merge(temp_dat, pathways, by.x = "Hugo_Symbol", by.y = "Gene", all.x = TRUE)

  temp_dat = temp_dat[order(-AlteredSamples, -Pathway, Hugo_Symbol)]
  temp_dat = rbind(temp_dat[!Pathway %in% "Unknown"], temp_dat[Pathway %in% "Unknown"])

  nsamps = as.numeric(m@summary[ID == 'Samples', summary])

  temp_dat[, pctAltered := paste0(round(AlteredSamples/nsamps * 100, 2), "%")]
  temp_dat = temp_dat[,.(Hugo_Symbol, AlteredSamples, pctAltered, Pathway, tcga_driver)]
  temp_dat
}
