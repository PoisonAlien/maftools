#' extract gene summary from MAF or GISTIC object
#' @name getGeneSummary
#' @rdname getGeneSummary
#' @param x An object of class MAF or GISTIC
#' @return gene summary table
#' @exportMethod getGeneSummary
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' getGeneSummary(laml)
setGeneric(name = "getGeneSummary", function(x) standardGeneric("getGeneSummary"))

#' extract sample summary from MAF or GISTIC object
#' @name getSampleSummary
#' @rdname getSampleSummary
#' @param x An object of class MAF or GISTIC
#' @return sample summary table
#' @exportMethod getSampleSummary
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' getSampleSummary(x = laml)
setGeneric(name = "getSampleSummary", function(x) standardGeneric("getSampleSummary"))

#' extract cytoband summary from GISTIC object
#' @name getCytobandSummary
#' @rdname getCytobandSummary
#' @param x An object of class GISTIC
#' @return summarizied gistic results by altered cytobands.
#' @exportMethod getCytobandSummary
#' @examples
#' all.lesions <- system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
#' amp.genes <- system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
#' del.genes <- system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
#' scores.gistic <- system.file("extdata", "scores.gistic", package = "maftools")
#' laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gistic)
#' getCytobandSummary(laml.gistic)
setGeneric(name = "getCytobandSummary", function(x) standardGeneric("getCytobandSummary"))

#' extract available fields from MAF object
#' @name getFields
#' @rdname getFields
#' @param x An object of class MAF
#' @return Field names in MAF file
#' @exportMethod getFields
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' getFields(x = laml)
setGeneric(name = "getFields", function(x) standardGeneric("getFields"))

#' extract annotations from MAF object
#' @name getAnnotations
#' @rdname getAnnotations
#' @param x An object of class MAF
#' @return annotations associated with samples in MAF
#' @exportMethod getAnnotations
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' getAnnotations(x = laml)
setGeneric(name = "getAnnotations", function(x) standardGeneric("getAnnotations"))

## Accessor methods
#' @rdname getGeneSummary
#' @aliases getGeneSummary
setMethod(f = "getGeneSummary",signature = "MAF", function(x) x@gene.summary)

## Accessor methods
#' @rdname getAnnotations
#' @aliases getAnnotations
setMethod(f = "getAnnotations",signature = "MAF", function(x) x@sample.anno)

#' @rdname getGeneSummary
#' @aliases getGeneSummary
setMethod(f = "getGeneSummary",signature = "GISTIC", function(x) x@gene.summary)

#' @rdname getSampleSummary
#' @aliases getSampleSummary
setMethod(f = "getSampleSummary",signature = "MAF", function(x) x@variant.classification.summary)

#' @rdname getSampleSummary
#' @aliases getSampleSummary
setMethod(f = "getSampleSummary",signature = "GISTIC", function(x) x@cnv.summary)

#' @rdname getCytobandSummary
#' @aliases getCytobandSummary
setMethod(f = "getCytobandSummary",signature = "GISTIC", function(x) x@cytoband.summary)

#' @rdname getFields
#' @aliases getFields
setMethod(f = "getFields",signature = "MAF", function(x) colnames(x@data))

