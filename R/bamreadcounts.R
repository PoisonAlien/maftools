#' Extract nucleotide counts for targeted variants from BAM files
#' @description bamreadcounts queries a BAM file at specific genomic positions and returns nucleotide counts (A, T, G, C) as well as insertion and deletion counts at each position. Function name is an homage to https://github.com/genome/bam-readcount
#' @param bam Path to indexed BAM file. BAM index (.bai) must exist in the same directory.
#' @param positions Either a 2-column data.frame with columns 'chr' and 'pos', or path to a TSV file with the same format.
#' @param mapq Minimum mapping quality threshold for reads. Default 10.
#' @param flag SAM flags to filter out. Default 1024 (duplicates).
#' @param fasta Path to indexed reference fasta file. If provided, output will include reference base at each position. Default NULL.
#' @param include_idxstats Calculate total mapped reads from BAM index. Default FALSE. When TRUE, result includes 'total_mapped_reads' attribute.
#' @param verbose Print progress messages. Default TRUE.
#' @return A data.frame with columns: chr, pos, A, T, G, C, INS, DEL, total_reads, coverage. If fasta is provided, includes additional ref_base column. If include_idxstats is TRUE, includes 'total_mapped_reads' attribute.
#' @examples
#' \dontrun{
#' # Using a data.frame
#' positions <- data.frame(chr = c("chr1", "chr2"), pos = c(12345, 67890))
#' result <- bamreadcounts("sample.bam", positions)
#' 
#' # Using a TSV file
#' result <- bamreadcounts("sample.bam", "positions.tsv")
#' 
#' # With reference fasta
#' result <- bamreadcounts("sample.bam", positions, fasta = "reference.fa")
#' 
#' # With idxstats for total mapped reads
#' result <- bamreadcounts("sample.bam", positions, include_idxstats = TRUE)
#' total_mapped <- attr(result, "total_mapped_reads")
#' }
#' @export
bamreadcounts <- function(bam, positions, mapq = 10, flag = 1024, fasta = NULL, include_idxstats = FALSE, verbose = TRUE) {
  
  # Validate inputs
  if (verbose) cat("Validating inputs...\n")
  .validate_bamreadcounts_inputs(bam, positions, mapq, flag)
  
  # Process positions input
  if (verbose) cat("Processing positions...\n")
  pos_df <- .process_positions_input(positions, verbose)
  
  # Call C function
  result <- .Call("bamrc_c", bam, pos_df$chr, pos_df$pos, as.integer(mapq), as.integer(flag), fasta, include_idxstats, verbose, NULL)
  
  if (verbose) cat("Done!\n")
  return(result)
}