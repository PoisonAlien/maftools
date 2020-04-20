#' Set Operations for MAF objects
#' @param x the first `MAF` object.
#' @param y the second `MAF` object.
#' @param mafObj Return output as an `MAF` object. Default `TRUE`
#' @param refAltMatch Set operations are done by matching ref and alt alleles in addition to loci (Default). Id FALSE only loci (chr, start, end positions) are matched.
#' @param ... other parameters passing to `subsetMaf` for subsetting operations.
#' @rdname setMaf
#' @export
#' @return subset table or an object of class \code{\link{MAF-class}}. If no overlaps found returns `NULL`
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' x <- subsetMaf(maf = laml, tsb = c('TCGA-AB-3009'))
#' y <- subsetMaf(maf = laml, tsb = c('TCGA-AB-2933'))
#' setdiffMAF(x, y)
#' intersectMAF(x, y) #Should return NULL due to no common variants
setdiffMAF <- function(x, y, mafObj = TRUE, refAltMatch = TRUE, ...) {
 stopifnot(inherits(x, "MAF"), inherits(y, "MAF"))

  args = list(...)
  if(length(args) > 0) {
    x = subsetMaf(x, mafObj = FALSE, ...)
    y = subsetMaf(y, mafObj = FALSE, ...)
  }

  if(refAltMatch){
    maf_x = data.table::rbindlist(l = list(nonsyn = x@data, syn = x@maf.silent), use.names = TRUE, fill = TRUE, idcol = "maf_slot")
    data.table::setkey(x = maf_x, Chromosome, Start_Position, End_Position)
    maf_y = data.table::rbindlist(l = list(y@data, y@maf.silent), use.names = TRUE, fill = TRUE)[,.(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2)]
    data.table::setkey(x = maf_y, Chromosome, Start_Position, End_Position)
    maf_x[, variant_ID := paste(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ":")]
    maf_y[, variant_ID := paste(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ":")]
  }else{
    maf_x = data.table::rbindlist(l = list(nonsyn = x@data, syn = x@maf.silent), use.names = TRUE, fill = TRUE, idcol = "maf_slot")
    data.table::setkey(x = maf_x, Chromosome, Start_Position, End_Position)
    maf_y = data.table::rbindlist(l = list(y@data, y@maf.silent), use.names = TRUE, fill = TRUE)[,.(Chromosome, Start_Position, End_Position)]
    data.table::setkey(x = maf_y, Chromosome, Start_Position, End_Position)
    maf_x[, variant_ID := paste(Chromosome, Start_Position, End_Position, sep = ":")]
    maf_y[, variant_ID := paste(Chromosome, Start_Position, End_Position, sep = ":")]
  }

  #Use faster character in vector operation
  maf_x_unique = maf_x[!maf_x$variant_ID %chin% maf_y$variant_ID]

  if (nrow(maf_x_unique) == 0) {
    warning("No X specific entries found!")
    return(NULL)
  }

  maf_x_unique[,variant_ID := NULL]
  maf_x_unique = droplevels.data.frame(maf_x_unique)

  if (!mafObj) {
    maf_x_unique
  } else {
    maf_x_unique = split(maf_x_unique, f = maf_x_unique$maf_slot)
    maf_x_unique[['syn']][,maf_slot := NULL]
    maf_x_unique[['nonsyn']][,maf_slot := NULL]
    mafSummary = summarizeMaf(maf = maf_x_unique[["nonsyn"]], anno = x@clinical.data, chatty = FALSE)

    maf_x_unique = MAF(data = maf_x_unique[['nonsyn']], variants.per.sample = mafSummary$variants.per.sample, variant.type.summary = mafSummary$variant.type.summary,
            variant.classification.summary = mafSummary$variant.classification.summary, gene.summary = mafSummary$gene.summary,
            summary = mafSummary$summary, maf.silent = maf_x_unique[['syn']], clinical.data = droplevels(mafSummary$sample.anno))
  }

  maf_x_unique
}


#' @rdname setMaf
#' @export
intersectMAF <- function(x, y, refAltMatch = TRUE, mafObj = TRUE, ...) {
  stopifnot(inherits(x, "MAF"), inherits(y, "MAF"))

  args = list(...)
  if(length(args) > 0) {
    x = subsetMaf(x, mafObj = FALSE, ...)
    y = subsetMaf(y, mafObj = FALSE, ...)
  }

  if(refAltMatch){
    maf_x = data.table::rbindlist(l = list(nonsyn = x@data, syn = x@maf.silent), use.names = TRUE, fill = TRUE, idcol = "maf_slot")
    data.table::setkey(x = maf_x, Chromosome, Start_Position, End_Position)
    maf_y = data.table::rbindlist(l = list(y@data, y@maf.silent), use.names = TRUE, fill = TRUE)[,.(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2)]
    data.table::setkey(x = maf_y, Chromosome, Start_Position, End_Position)
    maf_x[, variant_ID := paste(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ":")]
    maf_y[, variant_ID := paste(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ":")]
  }else{
    maf_x = data.table::rbindlist(l = list(nonsyn = x@data, syn = x@maf.silent), use.names = TRUE, fill = TRUE, idcol = "maf_slot")
    data.table::setkey(x = maf_x, Chromosome, Start_Position, End_Position)
    maf_y = data.table::rbindlist(l = list(y@data, y@maf.silent), use.names = TRUE, fill = TRUE)[,.(Chromosome, Start_Position, End_Position)]
    data.table::setkey(x = maf_y, Chromosome, Start_Position, End_Position)
    maf_x[, variant_ID := paste(Chromosome, Start_Position, End_Position, sep = ":")]
    maf_y[, variant_ID := paste(Chromosome, Start_Position, End_Position, sep = ":")]
  }

  #Use faster character in vector operation
  maf_x_common = maf_x[maf_x$variant_ID %chin% maf_y$variant_ID]

  if (nrow(maf_x_common) == 0) {
    warning("No common entries found!")
    return(NULL)
  }

  maf_x_common[,variant_ID := NULL]
  maf_x_common = droplevels.data.frame(maf_x_common)

  if (!mafObj) {
    maf_x_common
  } else {
    maf_x_common = split(maf_x_common, f = maf_x_common$maf_slot)
    maf_x_common[['syn']][,maf_slot := NULL]
    maf_x_common[['nonsyn']][,maf_slot := NULL]
    mafSummary = summarizeMaf(maf = maf_x_common[["nonsyn"]], anno = x@clinical.data, chatty = FALSE)

    maf_x_common = MAF(data = maf_x_common[['nonsyn']], variants.per.sample = mafSummary$variants.per.sample, variant.type.summary = mafSummary$variant.type.summary,
                       variant.classification.summary = mafSummary$variant.classification.summary, gene.summary = mafSummary$gene.summary,
                       summary = mafSummary$summary, maf.silent = maf_x_common[['syn']], clinical.data = droplevels(mafSummary$sample.anno))
  }

  maf_x_common

}

.mafSetKeys = function(maf){
  maf@data[,Chromosome := as.character(Chromosome)]
  maf@data[,Start_Position := as.numeric(as.character(Start_Position))]
  maf@data[,End_Position := as.numeric(as.character(End_Position))]

  maf@maf.silent[,Chromosome := as.character(Chromosome)]
  maf@maf.silent[,Start_Position := as.numeric(as.character(Start_Position))]
  maf@maf.silent[,End_Position := as.numeric(as.character(End_Position))]

  data.table::setkey(x = maf@data, Chromosome, Start_Position, End_Position)
  data.table::setkey(x = maf@maf.silent, Chromosome, Start_Position, End_Position)

  maf
}
