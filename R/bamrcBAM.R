.initialize_bam_access <- function(bam_path, verbose = TRUE) {

  # Load required library functions
  if (!requireNamespace("Rhtslib", quietly = TRUE)) {
    stop("Rhtslib package is required but not available. Please install it with: BiocManager::install('Rhtslib')")
  }

  tryCatch({
    # Open BAM file
    bam_file <- Rhtslib::BamFile(bam_path)

    # Test that we can read the header
    if (verbose) cat("Reading BAM header...\n")
    header_info <- Rhtslib::scanBamHeader(bam_file)

    # Extract chromosome information
    chr_info <- header_info[[1]]$targets
    if (length(chr_info) == 0) {
      stop("No chromosomes found in BAM header")
    }

    if (verbose) {
      cat(sprintf("BAM file contains %d chromosomes/contigs\n", length(chr_info)))
    }

    # Return connection object
    return(list(
      bam_file = bam_file,
      bam_path = bam_path,
      chr_info = chr_info,
      chr_names = names(chr_info)
    ))

  }, error = function(e) {
    stop("Failed to open BAM file: ", e$message)
  })
}

.process_bam_positions <- function(bam_connection, positions_df, mapq, flag, verbose) {

  # Initialize result data.frame
  result_df <- data.frame(
    chr = positions_df$chr,
    pos = positions_df$pos,
    A = integer(nrow(positions_df)),
    T = integer(nrow(positions_df)),
    G = integer(nrow(positions_df)),
    C = integer(nrow(positions_df)),
    INS = integer(nrow(positions_df)),
    DEL = integer(nrow(positions_df)),
    total_reads = integer(nrow(positions_df)),
    coverage = integer(nrow(positions_df)),
    stringsAsFactors = FALSE
  )

  # Process positions by chromosome for efficiency
  chr_groups <- split(seq_len(nrow(positions_df)), positions_df$chr)

  for (chr_name in names(chr_groups)) {
    chr_indices <- chr_groups[[chr_name]]
    chr_positions <- positions_df[chr_indices, ]

    if (verbose && length(chr_groups) > 1) {
      cat(sprintf("Processing chromosome %s (%d positions)...\n",
                  chr_name, length(chr_indices)))
    }

    # Check if chromosome exists in BAM
    if (!chr_name %in% bam_connection$chr_names) {
      if (verbose) {
        cat(sprintf("Warning: Chromosome '%s' not found in BAM file. Skipping %d positions.\n",
                    chr_name, length(chr_indices)))
      }
      next
    }

    # Process each position in this chromosome
    for (i in seq_along(chr_indices)) {
      row_idx <- chr_indices[i]
      pos <- chr_positions$pos[i]

      if (verbose && nrow(positions_df) > 100 && row_idx %% 100 == 0) {
        cat(sprintf("Processed %d/%d positions...\n", row_idx, nrow(positions_df)))
      }

      # Count nucleotides at this position
      counts <- .count_nucleotides_at_position(
        bam_connection, chr_name, pos, mapq, flag
      )

      # Store results
      result_df[row_idx, "A"] <- counts$A
      result_df[row_idx, "T"] <- counts$T
      result_df[row_idx, "G"] <- counts$G
      result_df[row_idx, "C"] <- counts$C
      result_df[row_idx, "INS"] <- counts$INS
      result_df[row_idx, "DEL"] <- counts$DEL
      result_df[row_idx, "total_reads"] <- counts$total_reads
      result_df[row_idx, "coverage"] <- counts$A + counts$T + counts$G + counts$C
    }
  }

  return(result_df)
}

.cleanup_bam_connection <- function(bam_connection) {
  tryCatch({
    if (!is.null(bam_connection$bam_file)) {
      close(bam_connection$bam_file)
    }
  }, error = function(e) {
    # Silently handle cleanup errors
  })

  invisible(NULL)
}
