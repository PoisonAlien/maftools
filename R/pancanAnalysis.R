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
#' @param pointSize size for scatter plot. Default 1.
#' @param labelSize label text size. Default 1
#' @return result table
#' @examples
#' laml.mutsig <- system.file("extdata", "LAML_sig_genes.txt.gz", package = "maftools")
#' pancanComparison(mutsigResults = laml.mutsig, qval = 0.1, cohortName = 'LAML', inputSampleSize = 200, label = 1)
#' @export

pancanComparison = function(mutsigResults, qval = 0.1, cohortName = 'input',
                            inputSampleSize = NULL, label = 1, genesToLabel = NULL, pointSize = 0.1, labelSize = 0.8){

  #Pancan results
  pancan = system.file('extdata', 'pancan.txt.gz', package = 'maftools')

  if(Sys.info()[['sysname']] == 'Windows'){
    pancan.gz = gzfile(description = pancan, open = 'r')
    pancan <- suppressWarnings( data.table::data.table(read.csv( file = pancan.gz, header = TRUE, sep = '\t', stringsAsFactors = FALSE)) )
    close(pancan.gz)
  } else{
    pancan = data.table::fread(cmd = paste('zcat <', pancan), sep = '\t', stringsAsFactors = FALSE)
  }

  #Input mutsig results
  if(as.logical(length(grep(pattern = 'gz$', x = mutsigResults, fixed = FALSE)))){
    #If system is Linux use fread, else use gz connection to read gz file.
    if(Sys.info()[['sysname']] == 'Windows'){
      mutsigResults.gz = gzfile(description = mutsigResults, open = 'r')
      suppressWarnings(mutsig <- data.table(read.csv(file = mutsigResults.gz, header = TRUE, sep = '\t', stringsAsFactors = FALSE, comment.char = '#')))
      close(mutsigResults.gz)
    } else{
      mutsig = suppressWarnings(data.table::fread(cmd = paste('zcat <', mutsigResults), sep = '\t', stringsAsFactors = FALSE, verbose = FALSE, data.table = TRUE, showProgress = TRUE, header = TRUE))
    }
  } else{
    mutsig = data.table::fread(input = mutsigResults, sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  }

  message(paste0('Significantly mutated genes in ', cohortName, ' (q < ', qval, '): ', nrow(mutsig[q < qval])))

  #0 q values (such as for TP53) are horrible for plotting since -log10(0) -> Inf. Changing them lowest machine epslion value for convenience
  mach.epsi = .Machine$double.eps
  mutsig$q = ifelse(test = as.numeric(as.character(mutsig$q)) == 0, yes = mach.epsi, no = mutsig$q)
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
  pancan.input.q[, log_q_pancan := -log10(pancan)]
  pancan.input.q[, log_q := -log10(q)]

  pancan.specific = pancan.input.q[q > qval & pancan < qval]
  pancan.rest = pancan.input.q[!gene %in% pancan.specific[,gene]]

  gTitle = paste0(cohortName, ' v/s Pan-cancer')
  xAxLab = paste0('-log10(', cohortName, ' q-value)')


  if(!is.null(genesToLabel)){
    lab_dat = pancan.rest[gene %in% genesToLabel]
  }else{
    if(label == 1){
      lab_dat = pancan.rest[q < qval & pancan > qval]
    }else if(label == 2){
      lab_dat = pancan.rest[q < qval & pancan < qval]
    }else if(label == 3){
      lab_dat = pancan.rest
    }else{
      stop('label can be 1, 2 or 3
         1: Labels genes mutated only in input cohort (Default)
         2: Labels genes mutated in both input and PanCancer cohort
         3: Labels all genes')
      }
    }

  message(paste0('Significantly mutated genes exclusive to ', cohortName, ' (q < ', qval, '): '))
  print(pancan.input.q[pancan > qval & q < qval])

  #return(pancan.input.q)
  par(mar = c(4, 4, 2, 2))
  bubble_plot(plot_dat = pancan.input.q, x_var = "log_q", y_var = "log_q_pancan", lab_dat = lab_dat,
              text_var = "gene", bubble_var = "log_q", text_size = labelSize)
  abline(v = -log10(qval), h = -log10(qval), col = 'maroon', lty = 2)
  title(paste0(cohortName, ' v/s Pan-cancer'), adj = 0)
  mtext(text =  paste0('-log10(', cohortName, ' q-value)'), side = 1, line = 2)
  mtext(text =  "-log10(Pan-can q-value)", side = 2, line = 2)

  return(pancan.input.q)

}
