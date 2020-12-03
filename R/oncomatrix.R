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
      return(NULL)
    }
  }

  if(add_missing){
    subMaf[, Hugo_Symbol := factor(x = Hugo_Symbol, levels = g)]
  }

  cnv_events = c(c("Amp", "Del"), as.character(subMaf[Variant_Type == "CNV"][, .N, Variant_Classification][, Variant_Classification]))
  cnv_events = unique(cnv_events)

  oncomat = data.table::dcast(data = subMaf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x, cnv = cnv_events){
                                #x = unique(as.character(x)) #>=2 distinct variant classification = Multi_Hit
                                x = as.character(x) # >= 2 same/distinct variant classification = Multi_Hit See #347
                                xad = x[x %in% cnv]
                                xvc = x[!x %in% cnv]

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

    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes, cnvc = cnv_events))
  }
}

update_vc_codes = function(om_op){
  uniq_vc = as.character(unique(unlist(as.numeric(unlist(apply(om_op$numericMatrix, 2, unique))))))
  missing_vc = uniq_vc[!uniq_vc %in% names(om_op$vc)]
  temp_names = names(om_op$vc)
  om_op$vc = c(om_op$vc,  rep("Complex_Event", length(missing_vc)))
  names(om_op$vc) = c(temp_names, rep(missing_vc, length(missing_vc)))
  om_op$vc
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

  if(is.numeric(anno[,1])){
    #anno.spl = split(anno, as.numeric(as.character(anno[,1]))) #sorting only first annotation
    mat_samps = colnames(numMat)[colnames(numMat) %in% rownames(anno)]
    anno = anno[mat_samps,, drop = FALSE]
    anno = anno[order(anno[,1], na.last = TRUE),, drop = FALSE]
    numMat.sorted = numMat[,rownames(anno)]
  }else{
    anno[,1] = ifelse(test = is.na(as.character(anno[,1])), yes = "NA", no = as.character(anno[,1])) #NAs are notorious; converting them to characters
    anno.spl = split(anno, as.factor(as.character(anno[,1]))) #sorting only first annotation

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
                       text_size = 1, col_var = NULL, return_dat = FALSE, showscale = FALSE, xlab = NULL, ylab = NULL){

  if(showscale){
    lo = layout(mat = matrix(data = c(1, 2), nrow = 1, ncol = 2), widths = c(5, 1))
  }
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
  axis(side = 2, at = y_ticks, las = 2, labels = abs(y_ticks))
  abline(h = y_ticks, v = x_ticks, lty = 2,
         col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.5), lwd = 0.75)

  mtext(text = xlab, side = 1, line = 2)
  mtext(text = ylab, side = 2, line = 3)

  if(!is.null(lab_dat)){
    text(x = lab_dat$x, y = lab_dat$y, labels = lab_dat$z_text, adj = 1, offset = 0.2, cex = text_size, col = lab_dat$color_var, xpd = TRUE)
  }

  if(showscale){
    par(mar = c(2, 1, 1, 1))
    plot(x = NA, ylim = c(0, 4.5), xlim = c(0, 2), axes = FALSE, xlab = NA, ylab = NA)
    bs = seq(min(plot_dat$size_z), max(plot_dat$size_z), length.out = 4)
    bslabs = round(seq(min(plot_dat$z), max(plot_dat$z), length.out = 4), 2)
    symbols(y = seq(3, 4, length.out = length(bs)), x = rep(0, length(bs)), circles = bs, add = TRUE, inches = 0.1, bg = "black", xpd = TRUE)
    text(x = rep(1, length(bs)), y = seq(3, 4, length.out = length(bs)), labels = bslabs, xpd = TRUE, adj =0)
    text(x = 0, y = 4.4, labels = "-log10(q)", xpd = TRUE, adj = 0)
  }
}


#Get plot layout for oncoplot
plot_layout = function(clinicalFeatures = NULL, drawRowBar = TRUE,
                       drawColBar = TRUE, draw_titv = FALSE, exprsTbl = NULL, legend_height = 4, anno_height = 1){

  if(is.null(clinicalFeatures)){
    if(draw_titv){
      if(!drawRowBar & !drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,rep(0,3)), nrow = 3, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 4, legend_height), widths = c(4,0.5))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 4, legend_height), widths = c(1, 4))
        }
      }else if(!drawRowBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3, 4,rep(0,4)), nrow = 4, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 12, 4, legend_height), widths = c(4, 0.5))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 12, 4, legend_height), widths = c(1, 4))
        }
      }else if(!drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 4, 4), widths = c(4, 1))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7,7), nrow = 3, ncol = 3, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 4, legend_height), widths = c(1, 4, 1))
        }
      }else{
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(4, 1), heights = c(4, 12, 4, legend_height))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,8,9,10,10,10), nrow = 4, ncol = 3, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(1,4, 1), heights = c(4, 12, 4, legend_height))
        }
      }
    }else{
      if(!drawRowBar & !drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,rep(0,2)), nrow = 2, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, legend_height), widths = c(4, 0.5))
        }else{
          mat_lo = matrix(data = c(1,2,3,3), nrow = 2, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, legend_height), widths = c(1, 4))
        }
      }else if(!drawRowBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,rep(0,3)), nrow = 3, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 12, legend_height), widths = c(4, 0.5))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 12, legend_height), widths = c(1, 4))
        }
      }else if(!drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,3), nrow = 2, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, legend_height), widths = c(4, 1))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,4,4), nrow = 2, ncol = 3, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, legend_height), widths = c(1, 4, 1))
        }
      }else{
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(4, 1), heights = c(4, 12, legend_height))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7,7), nrow = 3, ncol = 3, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(1, 4, 1), heights = c(4, 12, legend_height))
        }
      }
    }
  }else{
    if(draw_titv){
      if(!drawRowBar & !drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,rep(0,4)), nrow = 4, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, anno_height, 4, legend_height), widths = c(4, 0.5))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, anno_height, 4, legend_height), widths = c(1, 4))
        }
      }else if(!drawRowBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5,rep(0,5)), nrow = 5, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 12, anno_height, 4, legend_height), widths = c(4, 0.5))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,8,9,9), nrow = 5, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 12, anno_height, 4, legend_height), widths = c(1, 4))
        }
      }else if(!drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, anno_height, 4, legend_height), widths = c(4, 1))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,8,9,10,10,10), nrow = 4, ncol = 3, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, anno_height, 4, legend_height), widths = c(1, 4, 1))
        }
      }else{
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,8,9,9), nrow = 5, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(4, 1), heights = c(4, 12, anno_height, 4, legend_height))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,8,9,10,11,12,13,13,13), nrow = 5, ncol = 3, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(1, 4, 1), heights = c(4, 12, anno_height, 4, legend_height))
        }
      }
    }else{
      if(!drawRowBar & !drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,rep(0,3)), nrow = 3, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, anno_height, legend_height), widths = c(4, 0.5))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, anno_height, legend_height), widths = c(1, 4))
        }
      }else if(!drawRowBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,rep(0,4)), nrow = 4, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 12, anno_height, legend_height), widths = c(4, 0.5))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 12, anno_height, legend_height), widths = c(1, 4))
        }
      }else if(!drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, anno_height, legend_height), widths = c(4, 1))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7,7), nrow = 3, ncol = 3, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, anno_height, legend_height), widths = c(1, 4, 1))
        }
      }else{
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(4, 1), heights = c(4, 12, anno_height, legend_height))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,8,9,10,10,10), nrow = 4, ncol = 3, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(1, 4, 1), heights = c(4, 12, anno_height, legend_height))
        }
      }
    }
  }

  lo
}

get_domain_cols = function(){
  c("#f3a683", "#f7d794", "#778beb", "#e77f67", "#cf6a87", "#f19066",
    "#f5cd79", "#546de5", "#e15f41", "#c44569", "#786fa6", "#f8a5c2",
    "#63cdda", "#ea8685", "#596275", "#574b90", "#f78fb3", "#3dc1d3",
    "#e66767", "#303952")
}

get_vcColors = function(alpha = 1, websafe = FALSE, named = TRUE){
  if(websafe){
    col = c("#F44336", "#E91E63", "#9C27B0", "#673AB7", "#3F51B5", "#2196F3",
            "#03A9F4", "#00BCD4", "#009688", "#4CAF50", "#8BC34A", "#CDDC39",
            "#FFEB3B", "#FFC107", "#FF9800", "#FF5722", "#795548", "#9E9E9E",
            "#607D8B")
  }else{
    col = c(RColorBrewer::brewer.pal(11, name = "Paired"), RColorBrewer::brewer.pal(11,name = "Spectral")[1:3],'black', 'violet', 'royalblue', '#7b7060', '#535c68')
    col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  }

  if(named){
    names(col) = names = c('Nonstop_Mutation','Frame_Shift_Del','IGR','Missense_Mutation','Silent','Nonsense_Mutation',
                          'RNA','Splice_Site','Intron','Frame_Shift_Ins','In_Frame_Del','ITD','In_Frame_Ins',
                          'Translation_Start_Site',"Multi_Hit", 'Amp', 'Del', 'Complex_Event', 'pathway')
  }

  col
}

get_titvCol = function(alpha = 1){
  col = c("#F44336", "#3F51B5", "#2196F3", "#4CAF50", "#FFC107", "#FF9800")
  #col = c('coral4', 'lightcyan4', 'cornflowerblue', 'lightsalmon1', 'forestgreen', 'deeppink3')
  col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  names(col) = c('C>T', 'C>G', 'C>A', 'T>A', 'T>C', 'T>G')
  col
}

parse_annotation_dat = function(annotationDat = NULL, clinicalFeatures = NULL){

  if(class(annotationDat)[1] == 'MAF'){
    annotationDat = data.table::copy(getClinicalData(x = annotationDat))
    data.table::setDF(annotationDat)
  }else if(class(annotationDat)[1] %in% c('data.frame', 'data.table')){
    data.table::setDF(annotationDat)
  }else{
    return(NULL)
  }

  if(length(clinicalFeatures[!clinicalFeatures %in% colnames(annotationDat)]) > 0){
    message('Following columns are missing from annotation slot of MAF. Ignoring them..')
    print(clinicalFeatures[!clinicalFeatures %in% colnames(annotationDat)])
    clinicalFeatures = clinicalFeatures[clinicalFeatures %in% colnames(annotationDat)]
    if(length(clinicalFeatures) == 0){
      message('Make sure at-least one of the values from provided clinicalFeatures are present in annotation slot of MAF. Here are available annotaions..')
      print(colnames(annotationDat))
      stop('Zero annotaions to add! You can also provide custom annotations via annotationDat argument.')
    }
  }
  annotation = data.frame(row.names = annotationDat$Tumor_Sample_Barcode ,annotationDat[,clinicalFeatures, drop = FALSE], stringsAsFactors = FALSE)

  return(annotation)
}

get_anno_cols = function(ann, numericAnnoCol = NULL){
  ann_cols = list()
  if(is.null(numericAnnoCol)){
    numericAnnoCol =  RColorBrewer::brewer.pal(n = 9, name = "YlOrBr")
  }else{
    numericAnnoCol =  RColorBrewer::brewer.pal(n = 9, name = numericAnnoCol)
  }
  for(i in 1:ncol(ann)){
    if(is.numeric(ann[,i])){
      x = ann[,i]
      ann_lvls_cols = colorRampPalette(numericAnnoCol)(length(x))
      names(ann_lvls_cols) = x[order(x, na.last = TRUE)]
      ann_cols[[i]] = ann_lvls_cols
    }else{
      ann_lvls = unique(as.character(ann[,i]))
      if(length(ann_lvls) <= 9){
        ann_lvls_cols = RColorBrewer::brewer.pal(n = 9, name = 'Set1')[1:length(ann_lvls)]
      }else{
        ann_lvls_cols = colors()[sample(x = 1:100, size = length(ann_lvls), replace = FALSE)]
      }
      ann_cols[[i]] = ann_lvls_cols
      names(ann_cols[[i]]) = ann_lvls
    }
  }

  names(ann_cols) = colnames(ann)

  return(ann_cols)
}


pathway_load = function(maf){
  pathdb <- system.file("extdata", "BP_SMGs.txt.gz", package = "maftools")
  pathdb = data.table::fread(input = pathdb, skip = "Gene")
  pathdb = pathdb[!duplicated(Gene)][,.(Gene, Pathway)]
  pathdb$Pathway = gsub(pattern = " ", replacement = "_", x = pathdb$Pathway)
  pathdb_size = pathdb[,.N,Pathway]
  pathdb = split(pathdb, as.factor(pathdb$Pathway))


  altered_pws = lapply(pathdb, function(pw){
    x = suppressMessages(try(genesToBarcodes(maf = maf, genes = pw$Gene)))
    if(class(x) == "try-error"){
      pw_genes = NULL
    }else{
      pw_genes = names(genesToBarcodes(maf = maf, genes = pw$Gene, justNames = TRUE, verbose = FALSE))
    }
    pw_genes
  })

  mut_load = lapply(altered_pws, function(x){
    if(is.null(x)){
      nsamps =  0
    }else{
      nsamps = length(unique(as.character(unlist(
        genesToBarcodes(
          maf = maf,
          genes = x,
          justNames = TRUE
        )
      ))))
    }
    nsamps
  })

  altered_pws = as.data.frame(t(data.frame(lapply(altered_pws, length))))
  data.table::setDT(x = altered_pws, keep.rownames = TRUE)
  colnames(altered_pws) = c("Pathway", "n_affected_genes")
  altered_pws$Pathway = gsub(pattern = "\\.", replacement = "-", x = altered_pws$Pathway)
  altered_pws = merge(pathdb_size, altered_pws, all.x = TRUE)

  altered_pws[, fraction_affected := n_affected_genes/N]
  altered_pws$Mutated_samples = unlist(mut_load)
  nsamps = as.numeric(maf@summary[ID == "Samples", summary])
  altered_pws[,Fraction_mutated_samples := Mutated_samples/nsamps]
  altered_pws = altered_pws[order(Fraction_mutated_samples, fraction_affected, decreasing = TRUE)]

  altered_pws = altered_pws[!n_affected_genes %in% 0]
  altered_pws
}

#Update missing colors
update_colors = function(x, y){
  x = as.character(x)

  avail_colors = as.character(y[!names(y) %in% x])

  missing_entries = as.character(x)[!as.character(x) %in% names(y)]
  missing_entries = missing_entries[!missing_entries %in% ""]

  if(length(missing_entries) > 0){
    if(length(missing_entries) > length(avail_colors)){
      avail_colors = sample(x = colors(distinct = TRUE), size = length(missing_entries), replace = FALSE)
      names(avail_colors) = missing_entries
      y = c(y, avail_colors)
    }else{
      avail_colors = avail_colors[1:length(missing_entries)]
      names(avail_colors) = missing_entries
      y = c(y, avail_colors)
    }
  }
  y
}
