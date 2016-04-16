#' Plots segmented copy number data.
#' @details this function takes segmented copy number data and plots it. If MAF object is specified, all mutations are highlighted on the plot.
#' @param cbsFile CBS segmented copy number file. Column names should be Sample, Chromosome, Start, End, Num_Probes and Segment_Mean (log2 scale).
#' @param maf optional \code{\link{MAF}}
#' @param tsb If segmentation file contains many samples (as in gistic input), specify sample name here. Defualt plots all samples.
#' @param chr Just plot this chromosome.
#' @param savePlot If true plot is saved as pdf.
#' @param width width of plot
#' @param height height of plot
#' @param labelAll If true and if maf object is specified, maps all mutataions from maf onto segments. Default FALSE, maps only variants on copy number altered regions.
#' @param genes highlight only these variants
#' @return ggplot object
#' @export
#' @examples
#' tcga.ab.009.seg <- system.file("extdata", "TCGA.AB.009.hg19.seg.txt", package = "maftools")
#' plotCBSsegments(cbsFile = tcga.ab.009.seg)
#'

plotCBSsegments = function(cbsFile = NULL, maf = NULL, tsb = NULL, chr = NULL, savePlot = FALSE, width = 6, height = 3, labelAll = FALSE, genes = NULL){

  if(is.null(cbsFile)){
    stop('Required segmentation file !')
  }

  #Read segmentation file and change chromosome names
  seg = readSegs(seg = cbsFile)
  onlyContigs = as.character(seq(1:22))
  seg$Chromosome = gsub(pattern = 'chr', replacement = '', x = seg$Chromosome, fixed = TRUE)
  seg = seg[!Chromosome %in% c('X', 'Y')]
  seg = seg[Chromosome %in% onlyContigs]
  seg = seg[order(as.numeric(Chromosome))]
  setkey(x = seg, Chromosome, Start_Position, End_Position)

  #If user doesn't provide sample name
  if(is.null(tsb)){
    #Number of unique samples in segmentation file
    tsb = unique(seg[,Sample])
  }

  #If maf object is specified, map mutations on to segments
  if(!is.null(maf)){

    for(i in 1:length(tsb)){
      #Map mutations to segments
      tsb.mapped = mapMutsToSegs(seg = seg, maf = maf, tsb = tsb[i])

      if(!is.null(chr)){
        #If any specific chromosome is specified
        if(length(chr) >1){
          message('Multiple chromosomes specified. Using first entry.')
          chr = chr[1]
        }
        tsb.mapped = tsb.mapped[Chromosome == chr]
        p = plotCBSchr(segData = seg, tsb = tsb[i], chr = chr)

        #lable only variants on copy numbered regions
        if(!labelAll){

          if(!is.null(genes)){
            tsb.mapped.cn = tsb.mapped[Hugo_Symbol %in% genes]
          }else{
            tsb.mapped.cn = tsb.mapped[!CN >1.5 & CN < 2.5]
          }

          tsb.mapped.cn = tsb.mapped.cn[Chromosome == chr]
          p = p+geom_point(data = tsb.mapped, aes(x = Start_Position, y = Segment_Mean, label = Hugo_Symbol, color = 'green'), size = 0.6)+
            ggrepel::geom_text_repel(data = tsb.mapped.cn, aes(x = Start_Position, y = Segment_Mean,
            label = Hugo_Symbol), force = 10, nudge_y = 0.20, size = 2.5)+theme(legend.position = 'none')
        }else{
          #Lable all variants
          p = p+geom_point(data = tsb.mapped, aes(x = Start_Position, y = Segment_Mean, label = Hugo_Symbol, color = 'green'), size = 0.6)+
            ggrepel::geom_text_repel(data = tsb.mapped, aes(x = Start_Position, y = Segment_Mean,
            label = Hugo_Symbol), force = 10, nudge_y = 0.20, size = 2.5)+theme(legend.position = 'none')

        }
      }else{
        #Plot CBS segments for all chromosome
        p = suppressWarnings(plotCBS(segData = seg, tsb = tsb[i]))
        if(!labelAll){

          if(!is.null(genes)){
            tsb.mapped.cn = tsb.mapped[Hugo_Symbol %in% genes]
          }else{
            tsb.mapped.cn = tsb.mapped[!CN >1.5 & CN < 2.5]
          }

          p = p+geom_point(data = tsb.mapped, aes(x = Start_Position_updated, y = Segment_Mean, label = Hugo_Symbol, color = 'green'), size = 0.6)+
            ggrepel::geom_text_repel(data = tsb.mapped.cn, aes(x = Start_Position_updated, y = Segment_Mean,
            label = Hugo_Symbol), force = 10, nudge_y = 0.20, size = 2.5)+theme(legend.position = 'none')
        }else{
          p = p+geom_point(data = tsb.mapped, aes(x = Start_Position_updated, y = Segment_Mean, label = Hugo_Symbol), size = 0.6)+
            ggrepel::geom_text_repel(data = tsb.mapped, aes(x = Start_Position_updated, y = Segment_Mean,
            label = Hugo_Symbol), force = 10, nudge_y = 0.20, size = 2.5)+theme(legend.position = 'none')
        }
      }

      if(savePlot){
        cowplot::save_plot(filename = paste(tsb[i], '_segPlot.pdf', sep = ''), plot = p, base_height = height, base_width = width)
      }

    }

  }else{
    #Just plot segments
    for(i in 1:length(tsb)){
      if(!is.null(chr)){
        p = suppressWarnings(plotCBSchr(segData = seg, tsb = tsb[i], chr = chr))
      }else{
        p = suppressWarnings(plotCBS(segData = seg, tsb = tsb[i]))
      }

      if(savePlot){
        cowplot::save_plot(filename = paste(tsb[i], '_segPlot.pdf', sep = ''), plot = p, base_height = height, base_width = width)
      }
    }
  }

  return(p)
}
