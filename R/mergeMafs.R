#' Merge multiple maf files into a single MAF
#' @description Merges multiple maf files into a single MAF by matching column names.
#' @param mafs a vector of maf files.
#' @param MAFobj If TRUE, returns result as an \code{\link{MAF}} object. Default FALSE
#' @param ... additional arguments passed \code{\link{read.maf}}. Only applicable if MAFobj = TRUE.
#' @return data.table of merged MAFs or \code{\link{MAF}} object
#'
merge_mafs = function(mafs, MAFobj = FALSE, ...){

  maf = lapply(mafs, data.table::fread, stringsAsFactors = FALSE, fill = TRUE,
               showProgress = TRUE, header = TRUE, skip = "Hugo_Symbol")
  names(maf) = gsub(pattern = "\\.maf$", replacement = "", x = basename(path = mafs), ignore.case = TRUE)
  maf = data.table::rbindlist(l = maf, fill = TRUE, idcol = "sample_id", use.names = TRUE)

  if(MAFobj){
    maf = read.maf(maf = maf, ...)
  }

  maf
}
