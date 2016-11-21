#--- Small function to read segments
readSegs = function(seg){
  seg = suppressWarnings(fread(input = seg, sep = '\t', stringsAsFactors = FALSE, header = TRUE))
  seg$Sample = gsub(pattern = '-', replacement = '.', x = seg$Sample)

  #Replace chr x and y with numeric value (23 and 24) for better ordering
  seg$Chromosome = gsub(pattern = 'chr', replacement = '', x = seg$Chromosome, fixed = TRUE)
  seg$Chromosome = gsub(pattern = 'X', replacement = '23', x = seg$Chromosome, fixed = TRUE)
  seg$Chromosome = gsub(pattern = 'Y', replacement = '24', x = seg$Chromosome, fixed = TRUE)

  seg$Chromosome = as.character(seg$Chromosome)
  colnames(seg)[2:4] = c('Chromosome', 'Start_Position', 'End_Position')
  setkey(x = seg, Chromosome, Start_Position, End_Position)
  return(seg)
}

#--- Map mutaions from maf onto copynumber segments
mapMutsToSegs = function(seg, maf, tsb, build = 'hg19'){

  #onlyContigs = as.character(seq(1:22))

  seg.dat = seg[Sample %in% tsb]

  if(nrow(seg.dat) < 1){
    stop(paste('Sample',tsb, 'not found in segmentation file'))
  }

  seg.dat = transformSegments(segmentedData = seg.dat, build = build)
  setkey(x = seg.dat, Chromosome, Start_Position, End_Position)

  tsb.dat = subsetMaf(maf = maf, tsb = tsb, fields = 'Hugo_Symbol')

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

  tsb.dat = foverlaps(x = tsb.dat, y = seg.dat, by.x = c('Chromosome', 'Start_Position', 'End_Position'))
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
  p = ggplot(data = seg.dat)+geom_segment(data = seg.dat, aes(x = Start_Position, xend = End_Position, y = Segment_Mean, yend = Segment_Mean), color = 'maroon', size = 3)+
    cowplot::theme_cowplot(font_size = 8)+theme(axis.line.x = element_blank())+
    geom_hline(yintercept = 0, size = 0.5, linetype = 'dashed')+ggtitle(paste(tsb, ' : chr',chr, sep = ''))+xlab('Size')+ylab('Segment Mean')+
    cowplot::background_grid(major = 'onlyminor', minor = 'y')

  return(p)
}

#--- Plot complete segments
plotCBS = function(segData, tsb, build = 'hg19'){

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

  p = ggplot(data = seg.spl.transformed)+geom_segment(data = seg.spl.transformed, aes(x = Start_Position_updated, xend = End_Position_updated, y = Segment_Mean, yend = Segment_Mean, color = Chromosome), size = 3)+
    geom_vline(xintercept = chr.lens.sumsum[1:nchrs], linetype = 'dotted', size = 0.3)+
    cowplot::theme_cowplot(font_size = 8)+theme(legend.position = 'none')+xlab('Chromosome')+ylim(-2,2)+scale_x_continuous(breaks = chr.lens.sumsum[1:nchrs], labels = chr.labels[1:nchrs])+ylab('Segment Mean')+
    theme(axis.line.x = element_blank())+ggtitle(tsb)+cowplot::background_grid(major = 'onlyminor', minor = 'y')

  return(p)
}

# #--- Perform CBS segmentation on mutation data for raifall plot
# docbs = function(maf.snp){
#   x = DNAcopy::CNA(genomdat = maf.snp$diff, chrom = maf.snp$Chromosome, maploc = maf.snp$Start_Position, data.type = 'logratio', sampleid = tsb, presorted = TRUE)
#   CNA.smoothed <- DNAcopy::smooth.CNA(x)
#   segs = DNAcopy::segment(CNA.smoothed, verbose=0, min.width=2)
#   seg = DNAcopy::segs$output
#   colnames(seg) = c('Sample', 'Chromosome', 'Start_Position', 'End_Position', 'Num_Probes', 'Segment_Mean')
#   seg = transformSegments(segData = seg)
#   return(seg)
# }

