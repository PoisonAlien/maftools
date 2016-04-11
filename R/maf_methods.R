#' extract gene summary from MAF object
#' @name getGeneSummary
#' @rdname getGeneSummary
#' @param x An object of class MAF
#' @return gene summary table
#' @exportMethod getGeneSummary
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf, removeSilent = TRUE, useAll = FALSE)
#' getGeneSummary(laml)
setGeneric(name = "getGeneSummary", function(x) standardGeneric("getGeneSummary"))

#' extract sample summary from MAF object
#' @name getSampleSummary
#' @rdname getSampleSummary
#' @param x An object of class MAF
#' @return sample summary table
#' @exportMethod getSampleSummary
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf, removeSilent = TRUE, useAll = FALSE)
#' getSampleSummary(x = laml)
setGeneric(name = "getSampleSummary", function(x) standardGeneric("getSampleSummary"))

#' extract available fields from MAF object
#' @name getFields
#' @rdname getFields
#' @param x An object of class MAF
#' @return Field names in MAF file
#' @exportMethod getFields
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf, removeSilent = TRUE, useAll = FALSE)
#' getFields(x = laml)
setGeneric(name = "getFields", function(x) standardGeneric("getFields"))

## Accessor methods
#' @rdname getGeneSummary
#' @aliases getGeneSummary
setMethod(f = "getGeneSummary",signature = "MAF", function(x) x@gene.summary)

#' @rdname getSampleSummary
#' @aliases getSampleSummary
setMethod(f = "getSampleSummary",signature = "MAF", function(x) x@variant.classification.summary)

#' @rdname getFields
#' @aliases getFields
setMethod(f = "getFields",signature = "MAF", function(x) colnames(x@data))

