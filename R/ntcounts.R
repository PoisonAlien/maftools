#' extract nucleotide counts for targeted variants from the BAM file.
#' @description Given a BAM file and target loci, `ntcounts` fetches redcounts for A, T, G, C, Ins, and Del.
#' @param bam Input bam file. Required.
#' @param loci Loci file. First two columns should contain chromosome and position (1-based)
#' @param mapq Map quality. Default 10
#' @param sam_flag SAM FLAG to filter reads. Default 1024
#' @param op Output file basename. Default parses from BAM file
#' @param fa Indexed fasta file. If provided, extracts and adds reference base to the output tsv.
#' @useDynLib maftools, .registration = TRUE
#' @export

ntcounts = function(bam = NULL, loci = NULL, mapq = 10, sam_flag = 1024, op = NULL, fa = NULL){

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
    op_files = lapply(op_files, function(x) {
      paste0(x, ".tsv")
    })
  }else{
    if(length(op) != length(bam)){
      stop("No. of output file names must be equal to no. of BAM files.")
    }
    op_files = paste0(op, ".tsv")
  }

  if(is.null(fa)){
    fa = "NULL"
  }

  lapply(seq_along(bam), function(idx){
    withCallingHandlers(suppressWarnings(invisible(.Call("ntc", bam[idx], loci, mapq, sam_flag, fa, op[idx],  PACKAGE = "maftools"))))
  })

  res = lapply(seq_along(op_files), function(x){
    data.table::fread(file = op_files[[x]], sep = "\t", header = TRUE)
  })
  names(res) = op

  res
}
