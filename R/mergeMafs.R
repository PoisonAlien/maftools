#' Merge multiple mafs into single MAF
#' @description Merges multiple maf files/objects/data.frames into a single MAF.
#' @param mafs a list of \code{\link{MAF}} objects or data.frames or paths to MAF files.
#' @param verbose Default TRUE
#' @param ... additional arguments passed \code{\link{read.maf}}
#' @return \code{\link{MAF}} object
#'
merge_mafs = function(mafs, verbose = TRUE, ...){

  if(all(unlist(lapply(mafs, is, "MAF")))){

    if(verbose){
      cat(paste0("Merging ", length(mafs) ," MAF objects\n"))
    }
    mafs_dat = lapply(mafs, function(m) {
      data.table::rbindlist(l = list(m@data, m@maf.silent), use.names = TRUE, fill = TRUE)
    })
    mafs_dat = data.table::rbindlist(l = mafs_dat, fill = TRUE, use.names = TRUE)

    mafs_clin = lapply(mafs, function(m) {
      m@clinical.data
    })
    mafs_clin = data.table::rbindlist(l = mafs_clin, fill = TRUE, use.names = TRUE)
    maf = read.maf(maf = mafs_dat, clinicalData = mafs_clin, verbose = verbose, ...)
  }else if(all(unlist(lapply(mafs, is, "data.frame")))){
    if(verbose){
      cat(paste0("Merging ", length(mafs) ," MAF data.frames\n"))
    }
    mafs_dat = data.table::rbindlist(l = mafs, fill = TRUE, use.names = TRUE)
    maf = read.maf(maf = mafs_dat, verbose = verbose, ...)
  }else{
    if(verbose){
      cat(paste0("Merging ", length(mafs) ," MAF files\n"))
    }
    maf = lapply(mafs, data.table::fread, stringsAsFactors = FALSE, fill = TRUE,
                 showProgress = TRUE, header = TRUE, skip = "Hugo_Symbol")
    names(maf) = gsub(pattern = "\\.maf$", replacement = "", x = basename(path = unlist(mafs)), ignore.case = TRUE)
    maf = data.table::rbindlist(l = maf, fill = TRUE, idcol = "Source_MAF", use.names = TRUE)

    maf = read.maf(maf = maf, verbose = verbose, ...)
  }

  maf
}
