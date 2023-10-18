plotPathways = function(maf, pathways = NULL, fullPathway = FALSE, pathdb = NULL,
                                 removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue", fontSize = 0.6, showTumorSampleBarcodes = FALSE, sampleOrder = NULL, SampleNamefontSize = 0.6){


  if(is.null(pathdb)){
    pathdb <- system.file("extdata", "oncogenic_sig_patwhays.tsv", package = "maftools")
    pathdb = data.table::fread(file = pathdb)
  }else{
    if(is.data.frame(pathdb)){
      data.table::setDT(x = pathdb)
    }else{
      pathdb = data.table::fread(file = pathdb)
    }

    if(ncol(pathdb) >2){
      colnames(pathdb)[1:3] = c("Pathway", "Gene", "OG_TSG")
    }else{
      colnames(pathdb)[1:2] = c("Pathway", "Gene")
      pathdb$OG_TSG = ""
    }
  }

  #pathdb_size = pathdb[,.N,Pathway]
  pathdb = split(pathdb, as.factor(pathdb$Pathway))

  pathways = pathways[pathways %in% names(pathdb)]

  if(length(pathways) == 0){
    message("Available pathways..")
    print(names(pathdb))
    stop()
  }

  totSamps = as.numeric(maf@summary[3,summary])
  tsbs = levels(getSampleSummary(x = maf)[,Tumor_Sample_Barcode])

  par(mfrow = c(length(pathways), 1), bty="n",
      mar=c(2,4,2,2)+.1, las=1, tcl=-.25, cex=1)

  for(i in 1:length(pathways)){
    oncopath = pathdb[[pathways[i]]]
    data.table::setDF(x = oncopath, rownames = oncopath$Gene)
    oncopath$color_code = ifelse(test = oncopath$OG_TSG == "OG", yes = ogCol,
                                 no = ifelse(test = oncopath$OG_TSG == "TSG", yes = tsgCol, no = "black"))
    oncopath$color_code = ifelse(test = is.na(oncopath$color_code), yes = "black", no = oncopath$color_code)
    genes = oncopath$Gene
    path_mat = createOncoMatrix(m = maf, g = genes)
    if(is.null(path_mat)){
      stop("None of the genes in ", pathways[i], " pathway are mutated")
    }
    nm = path_mat$numericMatrix

    if(fullPathway){
      genes.missing = genes[!genes %in% rownames(nm)]
      if(length(genes.missing) > 0){
        genes.missing.numat = matrix(data = 0, ncol = ncol(nm), nrow = length(genes.missing))
        rownames(genes.missing.numat) = genes.missing
        nm = rbind(nm, genes.missing.numat)
      }
    }

    genes = rownames(nm)

    if(!removeNonMutated){
      tsb.include = matrix(data = 0, nrow = length(genes), ncol = length(tsbs[!tsbs %in% colnames(nm)]))
      colnames(tsb.include) = tsbs[!tsbs %in% colnames(nm)]
      rownames(tsb.include) = rownames(nm)
      nm = cbind(nm, tsb.include)
    }

    nm = t(apply(nm, 2, rev))
    nm[nm == 0] = NA
    nm[!is.na(nm)] = 1

    if(!is.null(sampleOrder)){
      sampleOrder = as.character(sampleOrder)
      sampleOrder = sampleOrder[sampleOrder %in% rownames(nm)]
      if(length(sampleOrder) == 0){
        stop("None of the provided samples are present in the input MAF")
      }
      nm = nm[sampleOrder, ,drop = FALSE]
    }

    if(showTumorSampleBarcodes){
      par(mar = c(4, 4, 2, 2), xpd = TRUE)
    }

    image(x = 1:nrow(nm), y = 1:ncol(nm), z = nm, axes = FALSE, xaxt="n", yaxt="n", xlab="", ylab="", col = "brown", ) #col = "#FC8D62"
    abline(h = (1:ncol(nm)) + 0.5, col = "white")
    abline(v = (1:nrow(nm)) + 0.5, col = "white")
    points(which(is.na(nm), arr.ind = TRUE), pch=".", col = "gray70")
    mtext(text = colnames(nm), side = 2, at = 1:ncol(nm), col = oncopath[colnames(nm), "color_code",],
          font = 3, line = 0.4, cex = fontSize)
    if(showTumorSampleBarcodes){
      # mtext(text = rownames(nm), side = 1, at = 1:nrow(nm),
      #       font = 3, line = 0.1, cex = SampleNamefontSize, las = 2)
      text(x =1:nrow(nm), y = par("usr")[3] - 0.2,
           labels = rownames(nm), srt = 45, font = 3, cex = SampleNamefontSize, adj = 1)
    }
    title(main = paste0(pathways[i], " pathway"), adj = 0)
  }
}
