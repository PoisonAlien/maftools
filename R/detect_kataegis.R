# Kategis detection algorithm by Moritz Goretzky @ WWU Munster, which exploits the definition of Kategis (six consecutive mutations with an avg. distance of 1000bp ) to idetify hyper mutated genomic loci.
# An algorithm which starts with a double-ended queue to which six consecutive mutations are added and their average intermutation distance is calculated.
# If the average intermutation distance is larger than 1000, one element is added at the back of the queue and one is removed from the front.
# If the average intermutation distance is less or equal to 1000, further mutations are added until the average intermutation distance is larger than 1000.
# After that all mutations in the double-ended queue are written into output as one kataegis and the double-ended queue is reinitialized with six mutations.

detect_kataegis_chr <- function(chr.dat) {
  chr.dat[, row_idx := 1:nrow(chr.dat)]

  start_idx <- 1
  end_idx <- 6
  queue <- chr.dat[start_idx:end_idx]

  kat_loc <- data.table::data.table()
  kat_id <- 1

  while (end_idx <= nrow(chr.dat)) {
    if (mean(diff(queue[, Start_Position], na.rm = TRUE), na.rm = TRUE) > 1000) {
      start_idx <- start_idx + 1
      end_idx <- end_idx + 1
      queue <- chr.dat[start_idx:end_idx]
    } else {
      while (mean(diff(queue[, Start_Position], na.rm = TRUE), na.rm = TRUE) <= 1000 &
        end_idx <= nrow(chr.dat)) {
        end_idx <- end_idx + 1
        queue <- chr.dat[start_idx:end_idx]
      }

      #---Summarize kat loci
      x <- chr.dat[(start_idx):c(end_idx - 1)] # start_idx not incremented after kat detected
      ycp <- data.table::data.table(
        Chromosome = unique(as.numeric(as.character(x[, Chromosome]))),
        Start_Position = x[, min(Start_Position)],
        End_Position = x[, max(Start_Position)],
        nMuts = nrow(x),
        Avg_intermutation_dist = mean(x[, diff(Start_Position)])
      )
      ycp[, Size := End_Position - Start_Position]

      ycp <- cbind(
        ycp,
        data.table::dcast(
          data = x[, .N, .(con.class, Tumor_Sample_Barcode)],
          Tumor_Sample_Barcode ~ con.class, value.var = "N"
        )
      )

      kat_loc <- data.table::rbindlist(l = list(kat_loc, ycp), fill = TRUE, use.names = TRUE)
      #---

      start_idx <- end_idx
      end_idx <- end_idx + 6
      queue <- chr.dat[(start_idx):(end_idx)]
    }
  }
  kat_loc
}


detect_kataegis <- function(maf.snp) {
  chr.spl <- split(maf.snp, maf.snp$Chromosome)
  kloci <- lapply(chr.spl, detect_kataegis_chr)
  kloci <- data.table::rbindlist(l = kloci, fill = TRUE, use.names = TRUE)
  if (nrow(kloci) > 0) {
    # kloci[,Filter := "PASS"]
    message("Kataegis detected at:")
    print(kloci)
    write.table(
      x = kloci,
      file = paste(unique(kloci[, Tumor_Sample_Barcode]), "Kataegis.tsv", sep = "_"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  kloci
}
