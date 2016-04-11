readSegs = function(seg){
  seg = suppressWarnings(fread(input = seg, sep = '\t', stringsAsFactors = FALSE, header = TRUE))
  seg$Sample = gsub(pattern = '-', replacement = '.', x = seg$Sample)
  seg = seg[order(seg$Chromosome)]
  seg$Chromosome = as.character(seg$Chromosome)
  colnames(seg)[2:4] = c('Chromosome', 'Start_Position', 'End_Position')
  setkey(x = seg, Chromosome, Start_Position, End_Position)
  return(seg)
}
