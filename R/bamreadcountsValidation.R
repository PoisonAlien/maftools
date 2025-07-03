.validate_bamreadcounts_inputs <- function(bam, positions, mapq, flag) {

  # Validate BAM file
  if (missing(bam) || is.null(bam) || !is.character(bam) || length(bam) != 1) {
    stop("bam must be a single character string specifying the BAM file path")
  }

  if (!file.exists(bam)) {
    stop("BAM file does not exist: ", bam)
  }

  # Check BAM file extension
  if (!grepl("\\.(bam|BAM)$", bam)) {
    stop("Input file must have .bam extension: ", bam)
  }

  # Check for BAM index
  bai_file <- paste0(bam, ".bai")
  if (!file.exists(bai_file)) {
    # Try alternative index naming
    bai_file_alt <- sub("\\.(bam|BAM)$", ".bai", bam)
    if (!file.exists(bai_file_alt)) {
      stop("BAM index file not found. Expected: ", bai_file, " or ", bai_file_alt)
    }
  }

  # Validate positions
  if (missing(positions) || is.null(positions)) {
    stop("positions must be provided as a data.frame or file path")
  }

  # Validate numeric parameters
  if (!is.numeric(mapq) || length(mapq) != 1 || mapq < 0) {
    stop("mapq must be a single non-negative number")
  }

  if (!is.numeric(flag) || length(flag) != 1 || flag < 0) {
    stop("flag must be a single non-negative integer")
  }

  invisible(TRUE)
}

.process_positions_input <- function(positions, verbose = TRUE) {

  if (is.character(positions) && length(positions) == 1) {
    # Input is a file path
    if (!file.exists(positions)) {
      stop("Positions file does not exist: ", positions)
    }

    if (verbose) cat("Reading positions from file:", positions, "\n")

    # Try to read the file
    tryCatch({
      if (requireNamespace("data.table", quietly = TRUE)) {
        pos_df <- data.table::fread(positions, header = TRUE, stringsAsFactors = FALSE)
        pos_df <- as.data.frame(pos_df)
      } else {
        # Fallback to base R
        pos_df <- read.table(positions, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      }
    }, error = function(e) {
      stop("Failed to read positions file: ", e$message)
    })

  } else if (is.data.frame(positions)) {
    # Input is already a data.frame
    pos_df <- as.data.frame(positions)

  } else {
    stop("positions must be either a data.frame or a file path")
  }

  # Validate data.frame structure
  if (nrow(pos_df) == 0) {
    stop("Positions data is empty")
  }

  # Check for required columns (case-insensitive)
  col_names <- tolower(names(pos_df))

  # Find chromosome column
  chr_col <- which(col_names %in% c("chr", "chromosome", "chrom", "seqnames"))
  if (length(chr_col) == 0) {
    stop("No chromosome column found. Expected column names: chr, chromosome, chrom, or seqnames")
  }
  if (length(chr_col) > 1) {
    chr_col <- chr_col[1]
    if (verbose) cat("Multiple chromosome columns found, using:", names(pos_df)[chr_col], "\n")
  }

  # Find position column
  pos_col <- which(col_names %in% c("pos", "position", "start", "coord"))
  if (length(pos_col) == 0) {
    stop("No position column found. Expected column names: pos, position, start, or coord")
  }
  if (length(pos_col) > 1) {
    pos_col <- pos_col[1]
    if (verbose) cat("Multiple position columns found, using:", names(pos_df)[pos_col], "\n")
  }

  # Create standardized data.frame
  result_df <- data.frame(
    chr = as.character(pos_df[[chr_col]]),
    pos = as.integer(pos_df[[pos_col]]),
    stringsAsFactors = FALSE
  )

  # Validate positions
  if (any(is.na(result_df$pos))) {
    stop("Some positions contain NA values")
  }

  if (any(result_df$pos < 1)) {
    stop("All positions must be positive integers (1-based coordinates)")
  }

  if (any(is.na(result_df$chr) | result_df$chr == "")) {
    stop("Some chromosome values are missing or empty")
  }

  # Remove duplicates and sort
  n_original <- nrow(result_df)
  result_df <- unique(result_df)
  if (nrow(result_df) < n_original && verbose) {
    cat(sprintf("Removed %d duplicate positions\n", n_original - nrow(result_df)))
  }

  # Sort by chromosome and position for efficiency
  result_df <- result_df[order(result_df$chr, result_df$pos), ]
  rownames(result_df) <- NULL

  if (verbose) cat(sprintf("Processed %d unique positions\n", nrow(result_df)))

  return(result_df)
}
