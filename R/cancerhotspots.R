#' Genotype known cancer hotspots from the tumor BAM file
#' @description `cancerhotspots` allows rapid genotyping of known somatic variants from the tumor BAM files. This facilitates to get a quick overlook of known somatic hot-spots in a matter of minutes, without spending hours on variant calling and annotation. In simple words, it fetches nucleotide frequencies of known somatic hotspots and prioritizes them based on allele frequency.
#' Output includes a browsable/sharable HTML report of candidate variants.
#' Known cancerhotspots for both GRCh37 and GRCh38 assemblies (3180 variants) are included. This should be sufficient and cover most of the known driver genes/events. See Reference for details.
#' @param bam Input bam file. Required.
#' @param refbuild Default "GRCh37". Can be "GRCh37", "GRCh38", "hg19", "hg38"
#' @param mapq Map quality. Default 10
#' @param sam_flag SAM FLAG to filter reads. Default 1024
#' @param vaf VAF threshold. Default 0.05 [Variant filter]
#' @param t_depth Depth of coverage threshold. Default 30 [Variant filter]
#' @param t_alt_count Min. number of reads supporting tumor allele . Default 8 [Variant filter]
#' @param op Output file basename. Default parses from BAM file
#' @param fa Indexed fasta file. If provided, extracts and adds reference base to the output tsv.
#' @param browse If TRUE opens the html file in browser
#' @useDynLib maftools, .registration = TRUE
#' @seealso \code{\link{cancerhotspotsAggr}}
#' @references Chang MT, Asthana S, Gao SP, et al. Identifying recurrent mutations in cancer reveals widespread lineage diversity and mutational specificity. Nat Biotechnol. 2016;34(2):155-163. doi:10.1038/nbt.3391
#' @export
cancerhotspots = function(bam = NULL, refbuild = "GRCh37", mapq = 10, sam_flag = 1024, vaf = 0.05, t_depth = 30, t_alt_count = 8,
               op = NULL, fa = NULL, browse = FALSE){

  if(any(is.null(bam), is.null(refbuild))){
    stop("Missing BAM or refbuild name!")
  }

  refbuild = match.arg(arg = refbuild, choices = c("GRCh37", "GRCh38", "hg19", "hg38"))

  bam_ext = substr(x = basename(bam), start = nchar(basename(path = bam))-3, nchar(basename(bam)))

  if(bam_ext != ".bam"){
    stop("Input file is not a BAM file: ", bam)
  }

  if(!file.exists(bam)){
    stop("BAM file does not exist: ", bam)
  }

  if(is.null(op)){
    op = gsub(pattern = "\\.bam$", replacement = "", x = basename(bam), ignore.case = TRUE)
  }

  if(is.null(fa)){
    fa = "NULL"
  }

  if(refbuild == "GRCh37"){
    bed = system.file("extdata", "cancerhotspots_v2_GRCh37.tsv", package = "maftools")
  }else if(refbuild == "GRCh38"){
    bed = system.file("extdata", "cancerhotspots_v2_GRCh38.tsv", package = "maftools")
  }else if(refbuild == "hg19"){
    bed = system.file("extdata", "cancerhotspots_v2_hg19.tsv", package = "maftools")
  }else{
    bed = system.file("extdata", "cancerhotspots_v2_hg38.tsv", package = "maftools")
  }


  withCallingHandlers(suppressWarnings(invisible(.Call("readb", bam, bed, t_depth, t_alt_count, vaf, fa, op, mapq, sam_flag, PACKAGE = "maftools"))))

  if(browse){
    browseURL(url = paste0(op, ".html"))
  }

  data.table::fread(file = paste0(op, ".tsv"))
}
