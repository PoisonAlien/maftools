createOncoMatrix = function(m, g = NULL, chatty = TRUE, add_missing = FALSE){

  if(is.null(g)){
    stop("Please provde atleast two genes!")
  }

  subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE, mafObj = FALSE)

  if(nrow(subMaf) == 0){
    if(add_missing){
      numericMatrix = matrix(data = 0, nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(numericMatrix) = g
      colnames(numericMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])

      oncoMatrix = matrix(data = "", nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(oncoMatrix) = g
      colnames(oncoMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])

      vc = c("")
      names(vc) = 0

      return(list(oncoMatrix = oncoMatrix, numericMatrix = numericMatrix, vc = vc))
    }else{

    }
    return(NULL)
  }

  if(add_missing){
    subMaf[, Hugo_Symbol := factor(x = Hugo_Symbol, levels = g)]
  }

  oncomat = data.table::dcast(data = subMaf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x){
                                x = unique(as.character(x))
                                xad = x[x %in% c('Amp', 'Del')]
                                xvc = x[!x %in% c('Amp', 'Del')]

                                if(length(xvc)>0){
                                  xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                }

                                x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                x = gsub(pattern = ';$', replacement = '', x = x)
                                x = gsub(pattern = '^;', replacement = '', x = x)
                                return(x)
                              } , value.var = 'Variant_Classification', fill = '', drop = FALSE)

  #convert to matrix
  data.table::setDF(oncomat)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[,-1, drop = FALSE])

  variant.classes = as.character(unique(subMaf[,Variant_Classification]))
  variant.classes = c('',variant.classes, 'Multi_Hit')
  names(variant.classes) = 0:(length(variant.classes)-1)

  #Complex variant classes will be assigned a single integer.
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)

  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes2)){
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }

  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }

  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)


  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    sampleId = colnames(mdf)
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId

    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId

    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]

    mdf = mdf[, -ncol(mdf)]

    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix

    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort

    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy

    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]

    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
}

#---- This is small function to sort genes according to total samples in which it is mutated.
sortByMutation = function(numMat, maf){

  geneOrder = getGeneSummary(x = maf)[order(MutatedSamples, decreasing = TRUE), Hugo_Symbol]
  numMat = numMat[as.character(geneOrder[geneOrder %in% rownames(numMat)]),, drop = FALSE]

  numMat[numMat != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
  tnumMat = t(numMat) #transposematrix
  numMat = t(tnumMat[do.call(order, c(as.list(as.data.frame(tnumMat)), decreasing = TRUE)), ]) #sort

  return(numMat)
}

#Thanks to Ryan Morin for the suggestion (https://github.com/rdmorin)
#original code has been changed with vectorized code, in-addition performs class-wise sorting.
sortByAnnotation <-function(numMat,maf, anno, annoOrder = NULL, group = TRUE, isNumeric = FALSE){
  anno[,1] = as.character(anno[,1])
  #anno[,1] = ifelse(test = is.na(anno[,1]), yes = "NA", no = anno[,1]) #NAs are notorious; converting them to characters
  #anno.spl = split(anno, anno[,1]) #sorting only first annotation (not converting to charcter)
  if(isNumeric){
    anno.spl = split(anno, as.numeric(as.character(anno[,1]))) #sorting only first annotation
  }else{
    anno[,1] = ifelse(test = is.na(anno[,1]), yes = "NA", no = anno[,1]) #NAs are notorious; converting them to characters
    anno.spl = split(anno, as.factor(as.character(anno[,1]))) #sorting only first annotation
  }


  anno.spl.sort = lapply(X = anno.spl, function(x){
    numMat[,colnames(numMat)[colnames(numMat) %in% rownames(x)], drop = FALSE]
  })

  if(group){
    #sort list according to number of elemnts in each classification
    anno.spl.sort = anno.spl.sort[names(sort(unlist(lapply(anno.spl.sort, ncol)), decreasing = TRUE, na.last = TRUE))]
  }

  if(!is.null(annoOrder)){
    annoSplOrder = names(anno.spl.sort)

    if(length(annoOrder[annoOrder %in% annoSplOrder]) == 0){
      message("Values in provided annotation order ", paste(annoOrder, collapse = ", ")," does not match values in clinical features. Here are the available features..")
      print(annoSplOrder)
      stop()
    }
    annoOrder = annoOrder[annoOrder %in% annoSplOrder]

    anno.spl.sort = anno.spl.sort[annoOrder]

    if(length(annoSplOrder[!annoSplOrder %in% annoOrder]) > 0){
      warning("Following levels are missing from the provided annotation order: ", paste(annoSplOrder[!annoSplOrder %in% annoOrder], collapse = ", "), immediate. = TRUE)
    }
  }

  numMat.sorted = c()
  for(i in 1:length(anno.spl.sort)){
    numMat.sorted  = cbind(numMat.sorted, anno.spl.sort[[i]])
  }

  return(numMat.sorted)
}

#Sort columns while keeping gene order
sortByGeneOrder = function(m, g){
  m = m[g,, drop = FALSE]
  tsbs= colnames(m)
  tsb.anno = data.table::data.table()

  for(i in 1:nrow(m)){
    x = m[i,]
    tsb.anno = rbind(tsb.anno, data.table::data.table(tsbs = names(x[x!=0]), gene = g[i]))
  }
  tsb.anno = tsb.anno[!duplicated(tsbs)]

  if(length(tsbs[!tsbs %in% tsb.anno[,tsbs]]) > 0){
    tsb.anno = rbind(tsb.anno, data.table::data.table(tsbs = tsbs[!tsbs %in% tsb.anno[,tsbs]], gene = NA))
  }
  data.table::setDF(tsb.anno)
  #m.sorted = sortByAnnotation(numMat = m, anno = tsb.anno)

  anno.spl = split(tsb.anno, as.factor(as.character(tsb.anno$gene)))
  anno.spl.sort = sapply(X = anno.spl, function(x){
    m[,colnames(m)[colnames(m) %in% x$tsbs], drop = FALSE]
  })

  numMat.sorted = c()
  for(i in 1:length(anno.spl.sort)){
    numMat.sorted  = cbind(numMat.sorted, anno.spl.sort[[i]])
  }

  return(numMat.sorted)
}

#plot_dat = plotting data
#lab_dat = data to be labelled
#x_var = x
#y_var = y
#bubble_var = z (variable name for bubble size)
#bubble_size (exclusive with bubble_var)
#text_size = font size for labels
#col_var = a vector color
bubble_plot = function(plot_dat, lab_dat = NULL, x_var = NULL, y_var = NULL,
                       bubble_var = NULL, bubble_size = 1, text_var = NULL,
                       text_size = 1, col_var = NULL, return_dat = FALSE){

  x_col_idx = which(colnames(plot_dat) == x_var)
  y_col_idx = which(colnames(plot_dat) == y_var)
  colnames(plot_dat)[c(x_col_idx, y_col_idx)] = c("x", "y")

  if(!is.null(col_var)){
    col_idx = which(colnames(plot_dat) == col_var)
    colnames(plot_dat)[col_idx] = c("color_var")
  }else{
    plot_dat$color_var = grDevices::adjustcolor("black", alpha.f = "0.75")
  }

  if(!is.null(lab_dat)){
    x_col_idx = which(colnames(lab_dat) == x_var)
    y_col_idx = which(colnames(lab_dat) == y_var)
    colnames(lab_dat)[c(x_col_idx, y_col_idx)] = c("x", "y")
    if(is.null(text_var)){
      stop("Missing text variable")
    }else{
      colnames(lab_dat)[which(colnames(lab_dat) == text_var)] = "z_text"
    }
    if(!is.null(col_var)){
      col_idx = which(colnames(lab_dat) == col_var)
      colnames(lab_dat)[col_idx] = c("color_var")
    }else{
      lab_dat$color_var = "black"
    }
  }

  if(!is.null(bubble_var)){

    if(bubble_var %in% c(x_var, y_var)){
      if(which(bubble_var == c(x_var, y_var)) == 2){
        plot_dat$size_z = sqrt(as.numeric(plot_dat$y)/pi)
      }else if(which(bubble_var == c(x_var, y_var)) == 1){
        plot_dat$size_z = sqrt(as.numeric(plot_dat$x)/pi)
      }
    }else{
      if(length(which(colnames(plot_dat) == bubble_var)) > 0){
        colnames(plot_dat)[which(colnames(plot_dat) == bubble_var)] = "z"
        plot_dat$size_z = sqrt(as.numeric(plot_dat$z)/pi)
      }else{
        plot_dat$size_z = bubble_size
      }
    }

    if(!is.null(lab_dat)){
      if(bubble_var %in% c(x_var, y_var)){
        if(which(bubble_var == c(x_var, y_var)) == 2){
          lab_dat$size_z = sqrt(as.numeric(lab_dat$y)/pi) * bubble_size
        }else if(which(bubble_var == c(x_var, y_var)) == 1){
          lab_dat$size_z = sqrt(as.numeric(lab_dat$x)/pi) * bubble_size
        }
      }else{
        colnames(lab_dat)[which(colnames(lab_dat) == bubble_var)] = "z"
        lab_dat$size_z = sqrt(as.numeric(lab_dat$z)/pi) * bubble_size
      }
    }
  }else{
    plot_dat$size_z = bubble_size
    if(!is.null(lab_dat)){
      lab_dat$size_z = bubble_size
    }
  }

  x_lims = as.integer(seq(min(as.numeric(plot_dat$x), na.rm = TRUE),
                          max(as.numeric(plot_dat$x), na.rm = TRUE), length.out = 4))
  x_lims[4] = as.integer(ceiling(max(as.numeric(plot_dat$x), na.rm = TRUE)))
  x_lims[1] = as.integer(floor(min(as.numeric(plot_dat$x), na.rm = TRUE)))
  x_ticks = pretty(x = x_lims, na.rm = TRUE)
  x_ticks[c(1, length(x_ticks))] = x_lims[c(1, 4)]

  y_lims = as.integer(seq(min(as.numeric(plot_dat$y), na.rm = TRUE),
                          max(as.numeric(plot_dat$y), na.rm = TRUE), length.out = 4))
  y_lims[4] = as.integer(ceiling(max(as.numeric(plot_dat$y), na.rm = TRUE)))
  y_lims[1] = as.integer(floor(min(as.numeric(plot_dat$y), na.rm = TRUE)))
  y_ticks = pretty(x = y_lims, na.rm = TRUE)
  y_ticks[c(1, length(y_ticks))] = y_lims[c(1, 4)]

  if(return_dat){
    return(list(plots_data = plot_dat, label_cords = lab_dat,
                x_lims = x_lims, y_lims = y_lims))
  }

  # plot(x = plot_dat$x, y = plot_dat$y, cex = plot_dat$size_z,
  #      pch = 16, col = plot_dat$color_var, axes = FALSE, xlim = x_lims[c(1, 4)],
  #      ylim = y_lims[c(1, 4)], xlab = NA, ylab = NA)
  suppressWarnings(symbols(x = plot_dat$x, y = plot_dat$y, circles = plot_dat$size_z, inches = 0.1, bg = plot_dat$color_var, xlim = x_lims[c(1, 4)],
          ylim = y_lims[c(1, 4)], xlab = NA, ylab = NA, fg = "white", axes = FALSE))
  axis(side = 1, at = x_ticks)
  axis(side = 2, at = y_ticks, las = 2)
  abline(h = y_ticks, v = x_ticks, lty = 2,
         col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.5), lwd = 0.75)

  if(!is.null(lab_dat)){
    # points(x = lab_dat$x, y = lab_dat$y, cex = lab_dat$size_z,
    #        pch = 16, col = lab_dat$color_var)
    symbols(x = lab_dat$x, y = lab_dat$y, circles = lab_dat$size_z,
            bg = lab_dat$color_var, add = TRUE, fg = "white", inches = 0.1)

    wordcloud::textplot(x = lab_dat$x, y = lab_dat$y, words = lab_dat$z_text,
                        cex = text_size, new = FALSE, show.lines = TRUE,
                        xlim = x_lims[c(1, 4)], ylim = y_lims[c(1, 4)], font = 3, col = lab_dat$color_var)
  }
}
