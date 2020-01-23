#--- Small function to read segments
readSegs = function(seg){
  seg = suppressWarnings(data.table::fread(input = seg, sep = '\t', stringsAsFactors = FALSE, header = TRUE))
  #seg$Sample = gsub(pattern = '-', replacement = '.', x = seg$Sample)

  #Replace chr x and y with numeric value (23 and 24) for better ordering
  seg$Chromosome = gsub(pattern = 'chr', replacement = '', x = seg$Chromosome, fixed = TRUE)
  seg$Chromosome = gsub(pattern = 'X', replacement = '23', x = seg$Chromosome, fixed = TRUE)
  seg$Chromosome = gsub(pattern = 'Y', replacement = '24', x = seg$Chromosome, fixed = TRUE)

  seg$Chromosome = as.character(seg$Chromosome)
  colnames(seg)[2:4] = c('Chromosome', 'Start_Position', 'End_Position')
  data.table::setkey(x = seg, Chromosome, Start_Position, End_Position)
  return(seg)
}

#--- Map mutations from maf onto copynumber segments
mapMutsToSegs = function(seg, maf, tsb, build = 'hg19'){

  #onlyContigs = as.character(seq(1:22))

  seg.dat = seg[Sample %in% tsb]

  if(nrow(seg.dat) < 1){
    stop(paste('Sample',tsb, 'not found in segmentation file'))
  }

  seg.dat = transformSegments(segmentedData = seg.dat, build = build)
  data.table::setkey(x = seg.dat, Chromosome, Start_Position, End_Position)

  tsb.dat = subsetMaf(maf = maf, tsb = tsb, fields = 'Hugo_Symbol', mafObj = FALSE)

  if(nrow(tsb.dat) < 1){
    stop(paste('Sample',tsb, 'not found in MAF'))
  }

  tsb.dat = tsb.dat[!Variant_Type %in% 'CNV']

  tsb.dat = tsb.dat[,.(Hugo_Symbol, Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode)]
  tsb.dat$Chromosome = gsub(pattern = 'X', replacement = '23', x = tsb.dat$Chromosome, fixed = TRUE)
  tsb.dat$Chromosome = gsub(pattern = 'Y', replacement = '24', x = tsb.dat$Chromosome, fixed = TRUE)
#   tsb.dat$Chromosome = gsub(pattern = 'chr', replacement = '', x = tsb.dat$Chromosome, fixed = TRUE)
  tsb.dat$Chromosome = as.character(tsb.dat$Chromosome)
#   tsb.dat = tsb.dat[!Chromosome %in% c('X', 'Y')]
#   tsb.dat = tsb.dat[Chromosome %in% onlyContigs]

  tsb.dat = data.table::foverlaps(x = tsb.dat, y = seg.dat, by.x = c('Chromosome', 'Start_Position', 'End_Position'))
  tsb.dat = tsb.dat[,.(Hugo_Symbol, Chromosome, i.Start_Position, i.End_Position,
                       Tumor_Sample_Barcode, Start_Position, End_Position, Segment_Mean, Start_Position_updated, End_Position_updated)]
  colnames(tsb.dat)[c(3:4, 6:7)] = c('Start_Position', 'End_Position', 'Segment_Start', 'Segment_End')

  suppressWarnings(tsb.dat[,CN := 2^(Segment_Mean)*2])
  return(tsb.dat)
}

#--- Change segment sizes into linear scale
transformSegments = function(segmentedData, build = 'hg19'){

  build.opts = c('hg19', 'hg18', 'hg38')

  if(!build %in% build.opts){
    stop('Available reference builds: hg18, hg19, hg38')
  }

  if(build == 'hg19'){
    chr.lens = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                 155270560, 59373566)
  } else if(build == 'hg18'){
    chr.lens = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
                 158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
                 114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
                 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
  } else if(build == 'hg38'){ #hg38
    chr.lens = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
  } else{
    stop('Available reference builds: hg18, hg19, hg38')
  }

  segmentedData[,Start_Position := as.numeric(as.character(Start_Position))]
  segmentedData[,End_Position := as.numeric(as.character(End_Position))]

  #Replace chr x and y with numeric value (23 and 24) for better ordering
  segmentedData$Chromosome = gsub(pattern = 'chr', replacement = '', x = segmentedData$Chromosome, fixed = TRUE)
  segmentedData$Chromosome = gsub(pattern = 'X', replacement = '23', x = segmentedData$Chromosome, fixed = TRUE)
  segmentedData$Chromosome = gsub(pattern = 'Y', replacement = '24', x = segmentedData$Chromosome, fixed = TRUE)

  segmentedData$Chromosome = factor(x = segmentedData$Chromosome, levels = 1:24, labels = 1:24)

  segmentedData = segmentedData[order(Chromosome, Start_Position, decreasing = FALSE)]

  seg.spl = split(segmentedData, segmentedData$Chromosome)

  seg.spl.transformed = seg.spl[[1]]
  if(nrow(seg.spl.transformed) > 0){
    seg.spl.transformed$Start_Position_updated = seg.spl.transformed$Start_Position
    seg.spl.transformed$End_Position_updated = seg.spl.transformed$End_Position
  }

  chr.lens.sumsum = cumsum(chr.lens)

  for(i in 2:length(seg.spl)){

    x.seg = seg.spl[[i]]
    if(nrow(x.seg) > 0){
      x.seg$Start_Position_updated = x.seg$Start_Position + chr.lens.sumsum[i-1]
      x.seg$End_Position_updated = x.seg$End_Position + chr.lens.sumsum[i-1]
    }
    seg.spl.transformed = rbind(seg.spl.transformed, x.seg, fill = TRUE)
  }

  return(seg.spl.transformed)
}

#--- Plot specific chromosome from segmentation data
plotCBSchr = function(segData, chr, tsb){
  seg = segData[Sample %in% tsb]

  if(nrow(seg) < 1){
    stop(paste('Sample',tsb, 'not found in segmentation file'))
  }

  seg.dat = seg[Chromosome == chr]

}

#--- Plot complete segments
plotCBS = function(segData, tsb, build = 'hg19', chr.colors = NULL, y_lims = NULL, rect_size = 0.1){

  segData = segData[Sample %in% tsb]

  if(nrow(segData) < 1){
    stop(paste('Sample',tsb, 'not found in segmentation file'))
  }

  segData = segData[order(Chromosome, Start_Position)]

  seg.spl.transformed = transformSegments(segmentedData = segData, build = build)

  if(build == 'hg19'){
    chr.lens = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                 155270560, 59373566)
  } else if(build == 'hg18'){
    chr.lens = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
                 158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
                 114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
                 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
  } else if(build == 'hg38'){
    chr.lens = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
  }

  chr.lens.sumsum = cumsum(chr.lens)
  nchrs = length(unique(seg.spl.transformed$Chromosome))

  chr.labels= c(1:22, 'X', 'Y')
  chr.lvls = levels(seg.spl.transformed[,Chromosome])

  print(chr.colors)
  if(is.null(chr.colors)){
    chr.colors = c('gray70', 'midnightblue')
  }

  chr.colors = rep(x = chr.colors, times = length(chr.lvls))[1:length(chr.lvls)]
  names(chr.colors) = chr.lvls

  if(is.null(y_lims)){
    y_lims = pretty(round(range(seg.spl.transformed$Segment_Mean, na.rm = TRUE, finite = TRUE), digits = 2))
  }

  par(mar = c(3, 3, 2, 1))
  plot(NA, xlim = c(0, max(seg.spl.transformed$End_Position_updated, na.rm = TRUE)),
       ylim = range(y_lims, na.rm = TRUE), axes = FALSE, ann = FALSE)
  axis(side = 1, at = chr.lens.sumsum, labels = chr.lvls, lwd.ticks = 1, tick = FALSE, line = -1)
  axis(side = 2,
       at = ,
       lwd.ticks = 1, las = 2)
  abline(v = chr.lens.sumsum, col = 'gray70', lwd = 0.25)
  rect(xleft = seg.spl.transformed$Start_Position_updated, xright =  seg.spl.transformed$End_Position_updated,
       ybottom = seg.spl.transformed$Segment_Mean-rect_size, ytop = seg.spl.transformed$Segment_Mean+rect_size,
       col = rep(chr.colors, seg.spl.transformed[,.N,Chromosome][,N]),
       border = rep(chr.colors, seg.spl.transformed[,.N,Chromosome][,N]))
  title(main = tsb, adj = 0, font.main = 3)
  mtext(text = "Chromosome", side = 1, line = 1, font = 3)
  mtext(text = "Segment mean", side = 2, line = 2, font = 3)
}
