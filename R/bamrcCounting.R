.count_nucleotides_at_position <- function(bam_connection, chr_name, position, mapq, flag) {

  # Initialize counts
  counts <- list(A = 0L, T = 0L, G = 0L, C = 0L, INS = 0L, DEL = 0L, total_reads = 0L)

  tryCatch({
    # Create genomic range for the position
    gr <- GenomicRanges::GRanges(
      seqnames = chr_name,
      ranges = IRanges::IRanges(start = position, end = position)
    )

    # Define what fields to extract from BAM
    what_fields <- c("seq", "cigar", "mapq", "flag", "pos")

    # Create scan parameters
    scan_param <- Rhtslib::ScanBamParam(
      which = gr,
      what = what_fields,
      flag = Rhtslib::scanBamFlag(
        isUnmappedQuery = FALSE,
        isSecondaryAlignment = FALSE
      )
    )

    # Scan BAM file
    bam_data <- Rhtslib::scanBam(bam_connection$bam_file, param = scan_param)[[1]]

    # Process each read
    if (length(bam_data$pos) > 0) {
      for (i in seq_along(bam_data$pos)) {

        # Apply quality filters
        read_mapq <- bam_data$mapq[i]
        read_flag <- bam_data$flag[i]

        if (is.na(read_mapq) || read_mapq < mapq) {
          next
        }

        if (!is.na(read_flag) && bitwAnd(read_flag, flag) != 0) {
          next
        }

        counts$total_reads <- counts$total_reads + 1L

        # Process this read
        read_counts <- .process_single_read(
          seq = bam_data$seq[[i]],
          cigar = bam_data$cigar[[i]],
          read_start = bam_data$pos[i],
          target_pos = position
        )

        # Add to total counts
        counts$A <- counts$A + read_counts$A
        counts$T <- counts$T + read_counts$T
        counts$G <- counts$G + read_counts$G
        counts$C <- counts$C + read_counts$C
        counts$INS <- counts$INS + read_counts$INS
        counts$DEL <- counts$DEL + read_counts$DEL
      }
    }

  }, error = function(e) {
    # Return zero counts on error
    warning("Error processing position ", chr_name, ":", position, " - ", e$message)
  })

  return(counts)
}

.process_single_read <- function(seq, cigar, read_start, target_pos) {

  # Initialize counts for this read
  read_counts <- list(A = 0L, T = 0L, G = 0L, C = 0L, INS = 0L, DEL = 0L)

  tryCatch({
    # Convert sequence to character if needed
    if (methods::is(seq, "DNAString")) {
      seq_char <- as.character(seq)
    } else {
      seq_char <- as.character(seq)
    }

    # Parse CIGAR string
    cigar_ops <- .parse_cigar_string(cigar)

    # Track positions
    ref_pos <- read_start  # Current reference position (1-based)
    read_pos <- 1         # Current position in read sequence (1-based)

    # Process each CIGAR operation
    for (op in cigar_ops) {
      op_type <- op$type
      op_length <- op$length

      if (op_type %in% c("M", "=", "X")) {
        # Match/mismatch operations

        # Check if target position falls within this operation
        if (target_pos >= ref_pos && target_pos < ref_pos + op_length) {
          # Target position is within this match block
          offset <- target_pos - ref_pos + 1  # 1-based offset within this block
          read_idx <- read_pos + offset - 1   # Index in read sequence

          if (read_idx <= nchar(seq_char)) {
            base <- substr(seq_char, read_idx, read_idx)
            base <- toupper(base)

            if (base == "A") read_counts$A <- 1L
            else if (base == "T") read_counts$T <- 1L
            else if (base == "G") read_counts$G <- 1L
            else if (base == "C") read_counts$C <- 1L
          }

          break  # Found the target position, exit loop
        }

        ref_pos <- ref_pos + op_length
        read_pos <- read_pos + op_length

      } else if (op_type == "I") {
        # Insertion

        # Check if insertion occurs exactly at target position
        if (ref_pos == target_pos) {
          read_counts$INS <- 1L
          break
        }

        read_pos <- read_pos + op_length
        # ref_pos stays the same for insertions

      } else if (op_type == "D") {
        # Deletion

        # Check if deletion spans the target position
        if (target_pos >= ref_pos && target_pos < ref_pos + op_length) {
          read_counts$DEL <- 1L
          break
        }

        ref_pos <- ref_pos + op_length
        # read_pos stays the same for deletions

      } else if (op_type == "N") {
        # Reference skip (intron for RNA-seq)
        ref_pos <- ref_pos + op_length

      } else if (op_type %in% c("S", "H")) {
        # Soft/hard clipping
        if (op_type == "S") {
          read_pos <- read_pos + op_length
        }
        # No change to ref_pos for clipping

      } else if (op_type == "P") {
        # Padding - no position changes
        next
      }

      # Early exit if we've passed the target position
      if (ref_pos > target_pos) {
        break
      }
    }

  }, error = function(e) {
    # Silently handle errors in individual read processing
  })

  return(read_counts)
}

.parse_cigar_string <- function(cigar) {

  if (is.na(cigar) || cigar == "" || cigar == "*") {
    return(list())
  }

  # Parse CIGAR using regular expression
  # CIGAR format: [0-9]+[MIDNSHPX=]
  pattern <- "([0-9]+)([MIDNSHPX=])"
  matches <- gregexpr(pattern, cigar, perl = TRUE)[[1]]

  if (matches[1] == -1) {
    return(list())
  }

  # Extract operations
  operations <- list()

  for (i in seq_along(matches)) {
    start_pos <- matches[i]
    match_length <- attr(matches, "match.length")[i]
    match_str <- substr(cigar, start_pos, start_pos + match_length - 1)

    # Extract length and type
    op_match <- regexpr(pattern, match_str, perl = TRUE)
    if (op_match == 1) {
      captures <- attr(op_match, "capture.start")
      capture_lengths <- attr(op_match, "capture.length")

      length_str <- substr(match_str, captures[1], captures[1] + capture_lengths[1] - 1)
      type_str <- substr(match_str, captures[2], captures[2] + capture_lengths[2] - 1)

      operations[[length(operations) + 1]] <- list(
        length = as.integer(length_str),
        type = type_str
      )
    }
  }

  return(operations)
}
