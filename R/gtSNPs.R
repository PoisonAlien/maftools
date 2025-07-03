#' Extract read counts from genetic markers for ASCAT analysis
#' @description The function will generate tsv files `<tumor/normal>_nucleotide_counts.tsv` that can be used for downstream analysis. Note that the function will process ~900K loci from Affymetrix Genome-Wide Human SNP 6.0 Array. The process can be sped up by increasing `nthreads` which will launch each chromosome on a separate thread. Currently hg19 and hg38 are supported. Files need to be further processed with \code{\link{prepAscat}} for tumor-normal pair, or \code{\link{prepAscat_t}} for tumor only samples.
#' @param t_bam Tumor BAM file. Required
#' @param n_bam Normal BAM file. Recommended
#' @param build Default hg19. Mutually exclusive with `loci`. Currently supported `hg19` and `hg38` and includes ca. 900K SNPs from Affymetrix Genome-Wide Human SNP 6.0 Array. SNP file has no `chr` prefix.
#' @param prefix Prefix to add or remove from contig names in loci file. For example, in case BAM files have `chr` prefix, set prefix = 'chr'
#' @param add If prefix is used, default is to add prefix to contig names in loci file. If false prefix will be removed from contig names.
#' @param mapq Minimum mapping quality. Default 10
#' @param sam_flag SAM FLAG to filter reads. Default 1024
#' @param loci A tab separated file with chr and position. If not available use `build` argument.
#' @param zerobased are coordinates zero-based. Default FALSE. Use only if `loci` is used.
#' @param op Output file basename. Default parses from BAM file
#' @param fa Indexed fasta file. If provided, extracts and adds reference base to the output tsv.
#' @param nthreads Number of threads to use. Default 4. Each chromosome will be launched on a separate thread. Works only on Unix and macOS.
#' @param verbose Default TRUE
#' @export
#' @seealso \code{\link{prepAscat}} \code{\link{prepAscat_t}} \code{\link{segmentLogR}}
#' @useDynLib maftools, .registration = TRUE
#' @import data.table

gtMarkers = function(t_bam = NULL, n_bam = NULL, build = "hg19", prefix = NULL, add = TRUE,
                      mapq = 10, sam_flag = 1024, loci = NULL, fa = NULL, op = NULL,
                      zerobased = FALSE, nthreads = 4, verbose = TRUE){

  if(is.null(t_bam)) stop("Missing tumor BAM file!")
  bam = c(t_bam)

  if(!is.null(n_bam)){
    bam = c(bam, n_bam)
  }


  if(is.null(loci)){
    if(build == "hg19"){
      download.file(url = "https://github.com/CompEpigen/ezASCAT/blob/main/inst/extdata/GRCh37_SNP6.tsv.gz?raw=true", destfile = "GRCh37_SNP6.tsv.gz")
      loci = "GRCh37_SNP6.tsv.gz"
    }else{
      download.file(url = "https://github.com/CompEpigen/ezASCAT/blob/main/inst/extdata/GRCh38_SNP6.tsv.gz?raw=true", destfile = "GRCh38_SNP6.tsv.gz")
      loci = "GRCh38_SNP6.tsv.gz"
    }
  }

  loci = data.table::fread(input = loci)
  colnames(loci)[1:2] = c("Chr", "start")

  if(!is.null(prefix)){
    if(add){
      loci$Chr = paste(prefix, loci$Chr, sep = '')
    }else{
      loci$Chr = gsub(pattern = prefix, replacement = '', x = loci$Chr, fixed = TRUE)
    }
  }

  data.table::setDF(x = loci)

  if(zerobased){
    loci$start = as.numeric(loci$start) + 1
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

  # Convert fa parameter for bamrc_c
  fasta_param = if(is.null(fa) || fa == "NULL") NULL else fa

  if(verbose){
    cat("Fetching readcounts from BAM files..\n")
  }

  res = list()
  bam_idxstats = list() #Store samtools idxstats

  for(b in bam){

    if(verbose){
      cat("Processing", basename(b), ":\n")
    }

    # Process all loci at once using bamrc_c with idxstats
    result = .Call("bamrc_c", b, loci$Chr, loci$start, as.integer(mapq), as.integer(sam_flag), fasta_param, TRUE, verbose, NULL, PACKAGE = "maftools")

    # Extract idxstats
    total_mapped = attr(result, "total_mapped_reads")
    idxstat = paste("#idxstats_mapped_reads", total_mapped)
    bam_idxstats[[length(bam_idxstats)+1]] = idxstat

    # Create loci column in format chr:pos
    result$loci = paste0(result$chr, ":", result$pos)

    # Reorder columns to match expected format
    if("ref_base" %in% colnames(result)) {
      result = result[, c("loci", "ref_base", "A", "T", "G", "C", "INS", "DEL")]
      colnames(result) = c("loci", "fa_ref", "A", "T", "G", "C", "Ins", "Del")
    } else {
      result = result[, c("loci", "A", "T", "G", "C", "INS", "DEL")]
      result$fa_ref = "N"  # Default if no fasta provided
      result = result[, c("loci", "fa_ref", "A", "T", "G", "C", "INS", "DEL")]
      colnames(result) = c("loci", "fa_ref", "A", "T", "G", "C", "Ins", "Del")
    }

    res[[length(res)+1]] = result
  }

  names(res) = op

  lapply(seq_along(res), function(idx){
    cat(paste0(bam_idxstats[[idx]], "\n"), file = paste0(op_files[[idx]], ".tsv"))
    data.table::fwrite(x = res[[idx]], file = paste0(op_files[[idx]], ".tsv"), append = TRUE, sep = "\t", na = "NA", quote = FALSE, col.names = TRUE)
    paste0(op_files[[idx]], ".tsv")
  }) |> unlist(use.names = FALSE)

}
