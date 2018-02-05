#' Perform PacCancer analysis
#' @description Takes MutSig results and compares them against PanCancer results.
#' @details This function takes MutSig results and compares them against panCancer cohort (~5000 tumor samples from 21 cancer types). This analysis can reveal
#' novel genes exclusively mutated in input cohort.
#'
#'@references Lawrence MS, Stojanov P, Mermel CH, et al. Discovery and saturation analysis of cancer genes across 21 tumor types. Nature. 2014;505(7484):495-501. doi:10.1038/nature12912.
#'
#' @param mutsigResults MutSig results (usually sig_genes.txt). Can be gz compressed.
#' @param qval qvalue threshold to define SMG. Default 0.1
#' @param cohortName Input cohort name.
#' @param inputSampleSize Sample size from MAF file used to generate mutSig results. Optional.
#' @param label Default 1. Can be 1, 2 or 3.
#' @param genesToLabel Default NULL. Exclusive with label argument.
#' @param normSampleSize normalizes gene sizes to draw bubble plot. Requires inputSampleSize. i.e, bubble sizes proportional to fraction of samples in which the gene is mutated.
#' @param pointSize size for scatter plot. Default 1.
#' @param labelSize label text size. Default 3
#' @param file basename for output file (both raw data and plot are saved)
#' @param width width of the file to be saved.
#' @param height height of the file to be saved.
#' @return ggplot object
#' @examples
#' laml.mutsig <- system.file("extdata", "LAML_sig_genes.txt.gz", package = "maftools")
#' pancanComparison(mutsigResults = laml.mutsig, qval = 0.1, cohortName = 'LAML', inputSampleSize = 200, label = 1, normSampleSize = TRUE)
#' @export

pancanComparison = function(mutsigResults, qval = 0.1, cohortName = 'input', inputSampleSize = NULL, label = 1, genesToLabel = NULL, normSampleSize = FALSE, file = NULL, width = 6, height = 6, pointSize = 3, labelSize = 3){

  #Pancan results
  pancan = system.file('extdata', 'pancan.txt.gz', package = 'maftools')

  if(Sys.info()[['sysname']] == 'Windows'){
    pancan.gz = gzfile(description = pancan, open = 'r')
    pancan <- suppressWarnings( data.table::data.table(read.csv( file = pancan.gz, header = TRUE, sep = '\t', stringsAsFactors = FALSE)) )
    close(pancan.gz)
  } else{
    pancan = data.table::fread(input = paste('zcat <', pancan), sep = '\t', stringsAsFactors = FALSE)
  }

  #Input mutsig results
  if(as.logical(length(grep(pattern = 'gz$', x = mutsigResults, fixed = FALSE)))){
    #If system is Linux use fread, else use gz connection to read gz file.
    if(Sys.info()[['sysname']] == 'Windows'){
      mutsigResults.gz = gzfile(description = mutsigResults, open = 'r')
      suppressWarnings(mutsig <- data.table(read.csv(file = mutsigResults.gz, header = TRUE, sep = '\t', stringsAsFactors = FALSE, comment.char = '#')))
      close(mutsigResults.gz)
    } else{
      mutsig = suppressWarnings(data.table::fread(input = paste('zcat <', mutsigResults), sep = '\t', stringsAsFactors = FALSE, verbose = FALSE, data.table = TRUE, showProgress = TRUE, header = TRUE))
    }
  } else{
    mutsig = data.table::fread(input = mutsigResults, sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  }

  message(paste0('Significantly mutated genes in ', cohortName, ' (q < ', qval, '): ', nrow(mutsig[q < qval])))

  #0 q values (such as for TP53) are horrible for plotting since -log10(0) -> Inf. Changing them lowest machine epslion value for convenience
  mach.epsi = .Machine$double.eps
  mutsig$q = ifelse(test = mutsig$q == 0, yes = mach.epsi, no = mutsig$q)
  mutsig.smg = mutsig[q < qval, gene]

  message(paste0('Significantly mutated genes in PanCan cohort (q <', qval, '): ', nrow(pancan[pancan < qval])))
  pancan.input.smg = unique(c(pancan[pancan < qval, gene], mutsig.smg))

  #message(paste0('Unique SMGs across ', cohortName ,' and PanCan: ', length(pancan.input.smg)))

  com.pancan.q = pancan[gene %in% pancan.input.smg, .(gene, pancan)]

  nNonCols = c('npat', 'nnon', 'n_nonsilent')

  if(length(nNonCols[nNonCols %in% colnames(mutsig)]) > 0){
    nNonCols = nNonCols[nNonCols %in% colnames(mutsig)][1]
    colnames(mutsig)[which(colnames(mutsig) == nNonCols)] = 'nMut'
  }else{
    stop('Column with number of non-silent mutations not found in input MutSigCV results!')
  }

  com.input.q = mutsig[gene %in% pancan.input.smg, .(gene, q, nMut)]

  pancan.input.q = merge(com.pancan.q, com.input.q, by = 'gene', all = TRUE)
  #Some genes are exlusively mutated in one cohort; such genes have q value set to NA. We will set them to 1 for plotting convenience
  pancan.input.q[is.na(pancan.input.q)] = 1

  pancan.specific = pancan.input.q[q > qval & pancan < qval]
  pancan.rest = pancan.input.q[!gene %in% pancan.specific[,gene]]

  #Plotting
  gTitle = paste0(cohortName, ' v/s Pan-cancer')
  xAxLab = paste0('-log10(', cohortName, ' q-value)')

  th = theme(legend.position = 'bottom', axis.text.x = element_text(face = "bold", size = 12),
             axis.title.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12),
             axis.title.y = element_text(face = "bold", size = 12),
             legend.title = element_text(face = "bold"), legend.text = element_text(size = 12, face = "bold"))


  if(normSampleSize){
    if(is.null(inputSampleSize)){
      stop('Missing inputSampleSize. Please provide number of samples from input MAF used to generate mutSig results.')
    }else{
      pancan.rest[,SampleFraction := nMut/inputSampleSize]
      pc.gg = ggplot(data = pancan.specific, aes(x = -log10(q), y = -log10(pancan), label = gene))+
        geom_point(alpha = 0.8, size = pointSize)+cowplot::theme_cowplot(font_size = 12, line_size = 1)+cowplot::background_grid(major = 'xy')+
        geom_hline(yintercept = -log10(qval), color = 'maroon', alpha = 0.8, linetype = 2)+
        geom_vline(xintercept = -log10(qval), color = 'maroon', alpha = 0.8, linetype = 2)+
        xlab(xAxLab)+ylab('-log10(Pan-can q-value)')+ggtitle(gTitle)+
        geom_point(data = pancan.rest, aes(x = -log10(q), y = -log10(pancan), size = SampleFraction), alpha  = 0.8)+th

    }
  }else{
    pc.gg = ggplot(data = pancan.specific, aes(x = -log10(q), y = -log10(pancan), label = gene))+
      geom_point(alpha = 0.8, size = pointSize)+cowplot::theme_cowplot(font_size = 12, line_size = 1)+cowplot::background_grid(major = 'xy')+
      geom_hline(yintercept = -log10(qval), color = 'maroon', alpha = 0.8, linetype = 2)+
      geom_vline(xintercept = -log10(qval), color = 'maroon', alpha = 0.8, linetype = 2)+
      xlab(xAxLab)+ylab('-log10(Pan-cancer q-value)')+ggtitle(gTitle)+
      geom_point(data = pancan.rest, aes(x = -log10(q), y = -log10(pancan)), size = pointSize)+
      th
  }

  message(paste0('Significantly mutated genes exclusive to ', cohortName, ' (q < ', qval, '): '))
  print(pancan.input.q[pancan > qval & q < qval])

  if(!is.null(genesToLabel)){
    pc.gg = pc.gg+ggrepel::geom_text_repel(data = pancan.rest[gene %in% genesToLabel], size = labelSize, color = 'red')
  }else{
    if(label == 1){
      pc.gg = pc.gg+ggrepel::geom_text_repel(data = pancan.rest[q < qval & pancan > qval], size = labelSize, color = 'red')
    }else if(label == 2){
      pc.gg = pc.gg+ggrepel::geom_text_repel(data = pancan.rest[q < qval & pancan < qval], size = labelSize, color = 'red')
    }else if(label == 3){
      pc.gg = pc.gg+ggrepel::geom_text_repel(data = pancan.rest, size = labelSize, color = 'red')
    }else{
      stop('label can be 1, 2 or 3
         1: Labels genes mutated only in input cohort (Default)
         2: Labels genes mutated in both input and PanCancer cohort
         3: Labels all genes')
    }

  }


  print(pc.gg)

  if(!is.null(file)){
    cowplot::save_plot(filename = paste0(file, '_pancan.pdf'), plot = pc.gg, base_height = height, base_width = width)
    pancan.input.q = pancan.input.q[order(q, pancan)]
    colnames(pancan.input.q)[2:3] = c('pancan_q', paste0(cohortName, '_q'))
    write.table(x = pancan.input.q, file = paste0(file, '_pancan.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
  }

  return(pc.gg)
}
