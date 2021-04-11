#' Display protein domains
#' @param gene HGNC symbol for which protein structure to be drawn.
#' @param refSeqID RefSeq transcript identifier for \code{gene} if known.
#' @param proteinID RefSeq protein identifier for \code{gene} if known.
#' @param domainAlpha Default 1
#' @param showLegend Default TRUE
#' @param bgBorderCol Default "black". Set to NA to remove.
#' @param axisTextSize text size x and y tick labels. Default c(1,1).
#' @param roundedRect Default TRUE. If `TRUE` domains are drawn with rounded corners. Requires \code{berryFunctions}
#' @param domainBorderCol Default "black". Set to NA to remove.
#' @param domainLabelSize text size for domain labels. Default 0.8
#' @param titleSize font size for title and subtitle. Default c(1.2, 1)
#' @param legendTxtSize Text size for legend. Default 0.8
#' @param legendNcol Default 1
#' @export
#' @examples
#' par(mfrow = c(2, 1))
#' plotProtein(gene = "KIT")
#' plotProtein(gene = "DNMT3A")
plotProtein = function(gene,
                       refSeqID = NULL,
                       proteinID = NULL,
                       domainAlpha = 0.9,
                       showLegend = FALSE,
                       bgBorderCol = "black",
                       axisTextSize = c(1, 1),
                       roundedRect = TRUE,
                       domainBorderCol = "black",
                       showDomainLabel = TRUE,
                       domainLabelSize = 0.8,
                       titleSize = c(1.2, 1),
                       legendTxtSize = 1,
                       legendNcol = 1
){
  if(is.null(gene)){
    stop('Please provide a gene name.')
  }

  geneID = gene
  #Protein domain source.
  gff = system.file('extdata', 'protein_domains.RDs', package = 'maftools')
  gff = readRDS(file = gff)
  data.table::setDT(x = gff)

  prot = gff[HGNC %in% geneID]

  if(nrow(prot) == 0){
    stop(paste('Structure for protein', geneID, 'not found.', sep=' '))
  }

  if(!is.null(refSeqID)){
    prot = prot[refseq.ID == refSeqID]
    if(nrow(prot) == 0){
      stop(paste0(refSeqID, " not found!"))
    }
  } else if(!is.null(proteinID)){
    prot = prot[protein.ID == proteinID]
    if(nrow(prot) == 0){
      stop(paste0(refSeqID, " not found!"))
    }
  } else{
    txs = unique(prot$refseq.ID)
    if(length(txs) > 1){
      message(paste(length(txs), ' transcripts available. Use arguments refSeqID or proteinID to manually specify tx name.', sep = ''))
      print(prot[!duplicated(protein.ID),.(HGNC, refseq.ID, protein.ID, aa.length)])
      prot = prot[which(prot$aa.length == max(prot$aa.length)),]
      if(length(unique(prot$refseq.ID)) > 1){
        prot = prot[which(prot$refseq.ID == unique(prot[,refseq.ID])[1]),]
        message(paste('Using longer transcript', unique(prot[,refseq.ID])[1], 'for now.', sep=' '))
      } else{
        message(paste('Using longer transcript', unique(prot[,refseq.ID])[1], 'for now.', sep=' '))
      }
    }
  }

  #Legth of protein
  len = as.numeric(max(prot$aa.length, na.rm = TRUE))
  #Remove NA's
  #prot = prot[!is.na(Label)]
  prot = prot[,domain_lenght := End - Start][order(domain_lenght, decreasing = TRUE)][,domain_lenght := NULL]

  xlimPos = pretty(0:max(prot$aa.length, na.rm = TRUE))
  xlimPos[length(xlimPos)] = max(prot$aa.length)

  domains = unique(prot[,Label])
  domain_cols = get_domain_cols()

  if(length(domains) > length(domain_cols)){
    domain_cols = sample(colours(), size = length(domains), replace = FALSE)
  }

  domain_cols = domain_cols[1:length(domains)]
  domain_cols = grDevices::adjustcolor(col = domain_cols, alpha.f = domainAlpha)
  names(domain_cols) = domains

  par(mar = c(0, 1, 2, 1))
  if(showDomainLabel){
    plot(0, 0, pch = NA, ylim = c(0.15, 0.85), xlim = c(0, len), axes = FALSE, xlab = NA, ylab = NA)
  }else{
    plot(0, 0, pch = NA, ylim = c(0, 0.85), xlim = c(0, len), axes = FALSE, xlab = NA, ylab = NA)
  }


  rect(xleft = 0, ybottom = 0.5, xright = len, ytop = 0.7, col = "#95a5a6", border = bgBorderCol)

  prot[, domainCol := domain_cols[prot[, Label]]]
  if(roundedRect){
    if(requireNamespace("berryFunctions", quietly = TRUE)){
      for(i in 1:nrow(prot)){
        berryFunctions::roundedRect(xleft = prot[i,Start], ybottom = 0.4, xright = prot[i,End], ytop = 0.8, col = prot[i, domainCol], border = domainBorderCol, rounding = 0.08)
      }
    }else{
      #warning("Package berryFunctions needed for roundedRect to work. Please install it and try again.")
      rect(xleft = prot[,Start], ybottom = 0.4, xright = prot[,End], ytop = 0.8, col = prot[,domainCol], border = domainBorderCol)
    }
  }else{
    rect(xleft = prot[,Start], ybottom = 0.4, xright = prot[,End], ytop = 0.8, col = prot[,domainCol], border = domainBorderCol)
  }

  text(x = xlimPos, y = 0.3, labels = xlimPos, col = "#34495e")
  #rect(xleft = xlimPos, ybottom = 0.22, xright = xlimPos, ytop = 0.24, col = "gray90")

  title(main = gene, adj = 0, font.main = 3, cex.main = titleSize[1], line = 0.8)
  title(main = unique(prot[,refseq.ID]), adj = 0, font.main = 1, line = -0.5, cex.main = titleSize[2])

  if(showDomainLabel){
    prot = prot[!duplicated(Label)]
    prot$pos = rowMeans(x = prot[,.(Start, End)], na.rm = FALSE)
    text(y = 0.6, x = prot$pos, labels = prot$Label, font = 3, cex = domainLabelSize)
  }else{
    legend(x = "bottomleft", legend = names(domain_cols), col = domain_cols, pch = 15, ncol = legendNcol, bty = "n", xpd = TRUE, cex = legendTxtSize)
  }
}
