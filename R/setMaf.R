#' Set Operations for MAF objects
#' @param x the first `MAF` object.
#' @param y the second `MAF` object.
#' @param ... other parameters passing to `subsetMaf` for subsetting operations.
#' @name setMaf
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' x <- subsetMaf(maf = laml, tsb = c('TCGA-AB-3009'))
#' y <- subsetMaf(maf = laml, tsb = c('TCGA-AB-2933'))
#' setdiffMAF(x, y)
#' intersectMAF(x, y)
#' # This is similar to merge_mafs()
#' unionMAF(x, y)
NULL


#' @rdname setMaf
setdiffMAF <- function(x, y, mafObj = TRUE, ...) {
 stopifnot(inherits(x, "MAF"), inherits(y, "MAF"))

  args = list(...)
  if (length(args) == 0) {
    maf_x = rbind(x@data, x@maf.silent)
    maf_y = rbind(y@data, y@maf.silent)
  } else {
    maf_x = subsetMaf(x, mafObj = FALSE, ...)
    maf_y = subsetMaf(y, mafObj = FALSE, ...)
  }

  maf_x[, Label := paste(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Type)]
  maf_y[, Label := paste(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Type)]

  set = setdiff(maf_x$Label, maf_y$Label)

  if (length(set) == 0) {
    message("No variants left.")
    return(NULL)
  }

  maf_x = maf_x[Label %in% set]
  maf_x$Label = NULL

  if (!mafObj) {
    maf_x
  } else {
    return(read.maf(maf_x, clinicalData = x@clinical.data, verbose = FALSE))
  }
}

#' @rdname setMaf
unionMAF <- function(x, y, mafObj = TRUE, ...) {
  stopifnot(inherits(x, "MAF"), inherits(y, "MAF"))

  args = list(...)

  if (length(args) == 0) {
    maf_x = rbind(x@data, x@maf.silent)
    maf_y = rbind(y@data, y@maf.silent)
  } else {
    maf_x = subsetMaf(x, mafObj = FALSE, ...)
    maf_y = subsetMaf(y, mafObj = FALSE, ...)
  }

  maf = unique(rbind(maf_x, maf_y, fill = TRUE))

  if (!mafObj) {
    return(maf)
  } else {
    cli_merged = rbind(x@clinical.data, y@clinical.data)
    return(read.maf(maf, clinicalData = cli_merged, verbose = FALSE))
  }

}


#' @rdname setMaf
#'
intersectMAF <- function(x, y, mafObj = TRUE, ...) {
  stopifnot(inherits(x, "MAF"), inherits(y, "MAF"))

  args = list(...)
  if (length(args) == 0) {
    maf_x = rbind(x@data, x@maf.silent)
    maf_y = rbind(y@data, y@maf.silent)
  } else {
    maf_x = subsetMaf(x, mafObj = FALSE, ...)
    maf_y = subsetMaf(y, mafObj = FALSE, ...)
  }

  maf_x[, Label := paste(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Type)]
  maf_y[, Label := paste(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Type)]

  set = intersect(maf_x$Label, maf_y$Label)
  if (length(set) == 0) {
    message("No common variants found.")
    return(NULL)
  }

  maf_x = maf_x[Label %in% set]
  maf_x$Label = NULL
  maf_y = maf_y[Label %in% set]
  maf_y$Label = NULL


  maf = rbind(maf_x, maf_y, fill = TRUE)
  if (!mafObj) {
    return(maf)
  } else {
    cli_merged = rbind(x@clinical.data, y@clinical.data)
    return(read.maf(maf, clinicalData = cli_merged, verbose = FALSE))
  }

}
