#' extract nucleotide counts for targeted variants from the BAM file.
#' @description Given a BAM file and target loci, `bamreadcounts` fetches redcounts for A, T, G, C, Ins, and Del. Function name is an homage to https://github.com/genome/bam-readcount
#' @param bam Input bam file(s). Required.
#' @param loci Loci file. Can be a tsv file or a data.frame. First two columns should contain chromosome and position (by default assumes coordinates are 1-based)
#' @param zerobased are coordinates zero-based. Default FALSE.
#' @param mapq Map quality. Default 10
#' @param sam_flag SAM FLAG to filter reads. Default 1024
#' @param op Output file basename. Default parses from BAM file
#' @param fa Indexed fasta file. If provided, extracts and adds reference base to the output tsv.
#' @param nthreads Number of threads to use. Each BAM file will be launched on a separate thread. Works only on Unix and macOS.
#' @useDynLib maftools, .registration = TRUE
#' @export

bamreadcounts = function(bam = NULL, loci = NULL, zerobased = FALSE, mapq = 10, sam_flag = 1024, op = NULL, fa = NULL, nthreads = 4){

  if(any(is.null(bam), is.null(loci))){
    stop("Missing BAM or loci file!")
  }

  op_files = lapply(bam, function(x){
    bam_ext = substr(x = basename(x), start = nchar(basename(path = x))-3, nchar(basename(x)))

    if(bam_ext != ".bam"){
      stop("Input file is not a BAM file: ", x)
    }

    if(!file.exists(x)){
      stop("BAM file does not exist: ", x)
    }
    gsub(pattern = "\\.bam$", replacement = "", x = basename(x), ignore.case = TRUE)
  })

  if(is.null(op)){
    op = as.character(unlist(op_files))
    op_files = lapply(op, function(x) {
      paste0(x, "_nucleotide_counts")
    })
    op_files = as.character(unlist(op_files))
  }else{
    if(length(op) != length(bam)){
      stop("No. of output file names must be equal to no. of BAM files.")
    }
    op_files = paste0(op, "_nucleotide_counts")
  }

  if(all(file.exists(op_files))){
    warning("Counts are already generated!")
    res = lapply(seq_along(op_files), function(x){
      data.table::fread(file = paste0(op_files[x], ".tsv"), sep = "\t", header = TRUE)
    })
    names(res) = op
    return(res)
  }

  if(is.null(fa)){
    fa = "NULL"
  }

  if(is.data.frame(loci)){
    colnames(loci)[c(1:2)] = c("chr", "start")
  }else if(file.exists(loci)){
    loci = data.table::fread(file = loci)
    colnames(loci)[c(1:2)] = c("chr", "start")
  }else{
    stop("loci must be a data.frame or a tsv file!")
  }
  data.table::setDF(x = loci)

  if(zerobased){
    loci$start = as.numeric(loci$start) + 1
  }

  lfile = tempfile(pattern = 'bamrc_', fileext = ".tsv")
  data.table::fwrite(x = loci[,c(1:2)], file = lfile, col.names = FALSE, sep = "\t", row.names = FALSE)

  cat("Fetching readcounts from BAM files..\n")
  parallel::mclapply(seq_along(bam), function(idx){
    cat("Processing", bam[idx], ":\n")
    withCallingHandlers(suppressWarnings(invisible(.Call("ntc", bam[idx], lfile, mapq, sam_flag, fa, op_files[idx],  PACKAGE = "maftools"))))
  }, mc.cores = nthreads)

  res = lapply(seq_along(op_files), function(x){
    data.table::fread(file = paste0(op_files[x], ".tsv"), sep = "\t", header = TRUE)
  })
  names(res) = op

  res
}
