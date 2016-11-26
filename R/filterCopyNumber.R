filterCopyNumber = function(seg, tsb.dat, tempCheck, tsb){
  if(nrow(seg) < 1){
    stop(paste('No copynumber data found for sample', tsb, sep = ' '))
  }else{

    #Overlap variants with segment data
    tsb.dat = data.table::foverlaps(x = tsb.dat, y = seg, by.x = c('Chromosome', 'Start_Position', 'End_Position'))
    tsb.dat = tsb.dat[,.(Hugo_Symbol, Chromosome, i.Start_Position, i.End_Position,
                         Tumor_Sample_Barcode, t_vaf, Start_Position, End_Position, Segment_Mean)]
    colnames(tsb.dat)[c(3:4, 7:8)] = c('Start_Position', 'End_Position', 'Segment_Start', 'Segment_End')

    #Convert log scale to absolute copynumber data
    suppressWarnings(tsb.dat[,CN := 2^(Segment_Mean)*2])

    if(nrow(tsb.dat[is.na(tsb.dat$CN)]) > 0){
      message(paste('Removed ', nrow(tsb.dat[is.na(tsb.dat$CN)]), ' variants with no copy number data.', sep = ''))
      print(tsb.dat[is.na(tsb.dat$CN)])
      tsb.dat = tsb.dat[!is.na(tsb.dat$CN)]
    }

    #Remove copy number altered variants.
    tsb.dat.cn.vars = tsb.dat[!CN >1.5 & CN < 2.5]
    if(nrow(tsb.dat.cn.vars) > 0){
      message('Copy number altered variants:')
      tsb.dat.cn.vars$cluster = 'CN_altered'
      print(tsb.dat.cn.vars)
      tempCheck = 1
    }
    tsb.dat = tsb.dat[CN >1.5 & CN < 2.5] #Copy number neutral variants
  }

  return(list(tsb.dat, tsb.dat.cn.vars, tempCheck))
}
