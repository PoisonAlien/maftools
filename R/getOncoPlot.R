getOncoPlot = function(maf, genes, removeNonMutated = FALSE, colors = NULL, showGenes = TRUE,
                       left = FALSE,
                       showTumorSampleBarcodes = FALSE, hmName = hmName, fs = 10, gfs = 10, tfs = 12,
                       clinicalFeatures = NULL, annotationColor = NULL, keepGeneOrder = FALSE,
                       includeSyn=FALSE){

  #-----preprocess matrix

  tsbs = levels(getSampleSummary(x = maf)[,Tumor_Sample_Barcode])

   om = createOncoMatrix(m = maf, g = genes)
  if(is.null(om)){
    nsamps = as.numeric(maf@summary[ID %in% 'Samples', summary])
    oncoMatrix = matrix(data = "", nrow = length(genes), ncol = nsamps)
    numericMatrix = matrix(data = 0, nrow = length(genes), ncol = nsamps)
    colnames(oncoMatrix) = as.character(getSampleSummary(maf)[,Tumor_Sample_Barcode])
    rownames(oncoMatrix) = as.character(genes)

    colnames(numericMatrix) = as.character(getSampleSummary(maf)[,Tumor_Sample_Barcode])
    rownames(numericMatrix) = as.character(genes)

    om = list(numericMatrix = numericMatrix, oncoMatrix =oncoMatrix)
  }



   mat_origin = om$numericMatrix
  numMat = om$numericMatrix


  if(ncol(mat_origin) < 2){
    stop('Cannot create oncoplot for single sample. Minimum two sample required ! ')
  }

  genes.missing = genes[!genes %in% rownames(mat_origin)]
  genes.present = genes[genes %in% rownames(mat_origin)]
  mat = mat_origin[genes.present,,drop = FALSE]

  if(nrow(mat) == 0){
    message(paste0('NOTE: Zero samples are mutated in ', hmName))
  }
  #remove nonmutated samples to improve visualization
  if(removeNonMutated){
    tsb = colnames(mat)
    tsb.exclude = colnames(mat[,colSums(mat) == 0, drop = FALSE])
    tsb.include = tsb[!tsb %in% tsb.exclude]
    mat = mat[,tsb.include, drop = FALSE]
  }

  if(keepGeneOrder){
    #Sort
    mat[mat != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmat = t(mat)
    mat_origin = om$oncoMatrix


     if(length(genes.missing) > 0){
        genes.missing.mat = t(matrix(data = '', ncol = ncol(mat), nrow = length(genes.missing)))
        genes.missing.numat = t(matrix(data = 0, ncol = ncol(mat), nrow = length(genes.missing)))
        colnames(genes.missing.mat) = genes.missing
        colnames(genes.missing.numat) = genes.missing
        mat = rbind(mat, t(genes.missing.mat))
        numMat = rbind(numMat, t(genes.missing.numat))
        #mat = mat[genes,]


        genes.missing.mat2 = t(matrix(data = '', ncol = ncol(mat_origin), nrow = length(genes.missing)))
        colnames(genes.missing.mat2) = genes.missing
        mat_origin = rbind(mat_origin, t(genes.missing.mat2))
     }
    char.mat = mat_origin[rownames(mat),, drop = FALSE]
    char.mat = char.mat[,colnames(mat), drop = FALSE]
    mat = char.mat
    numMat = numMat[rownames(mat),, drop = FALSE]
    numMat = numMat[,colnames(mat), drop = FALSE]
    numMat = sortByMutationKeepGeneOrder(numMat, genes)
  }else{
    #Sort
    mat[mat != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmat = t(mat)

    if(nrow(mat) == 1){
      g = rownames(mat)
      mat = t(tmat[do.call(order, c(as.list(as.data.frame(tmat)), decreasing = TRUE)), ,drop = FALSE])
      rownames(mat) = g
    }else{
      mat = t(tmat[do.call(order, c(as.list(as.data.frame(tmat)), decreasing = TRUE)), ])
    }

    mat_origin = om$oncoMatrix
    char.mat = mat_origin[rownames(mat),, drop = FALSE]
    char.mat = char.mat[,colnames(mat), drop = FALSE]
    mat = char.mat
    numMat = numMat[rownames(mat),, drop = FALSE]
    numMat = numMat[,colnames(mat), drop = FALSE]

    if(length(genes.missing) > 0){
      genes.missing.mat = t(matrix(data = '', ncol = ncol(mat), nrow = length(genes.missing)))
      genes.missing.numat = t(matrix(data = 0, ncol = ncol(mat), nrow = length(genes.missing)))
      colnames(genes.missing.mat) = genes.missing
      colnames(genes.missing.numat) = genes.missing
      mat = rbind(mat, t(genes.missing.mat))
      numMat = rbind(numMat, t(genes.missing.numat))
      #mat = mat[genes,]

      genes.missing.mat2 = t(matrix(data = '', ncol = ncol(mat_origin), nrow = length(genes.missing)))
      colnames(genes.missing.mat2) = genes.missing
      mat_origin = rbind(mat_origin, t(genes.missing.mat2))
    }
  }

  #final matrix for plotting
  mat = mat[genes , , drop = FALSE]
  mat_origin = mat_origin[genes , , drop = FALSE]
  numMat = numMat[genes,,drop = FALSE]


  #New version of complexheatmap complains about '' , replacing them with random strinf xxx

  mat[mat == ''] = 'xxx'

  #By default oncomatrix excludes non-mutated samples. Add rest here if user requests
  if(!removeNonMutated){
    tsb.include = matrix(data = '', nrow = length(genes), ncol = length(tsbs[!tsbs %in% colnames(mat)]))
    colnames(tsb.include) =tsbs[!tsbs %in% colnames(mat)]
    rownames(tsb.include) = rownames(mat)
    mat = cbind(mat, tsb.include)
  }

  #---------------------------------------Colors and vcs-------------------------------------------------

  if(is.null(colors)){
    col = c(RColorBrewer::brewer.pal(12,name = "Paired"), RColorBrewer::brewer.pal(11,name = "Spectral")[1:3],'black', 'red', 'green')
    names(col) = names = c('Nonstop_Mutation','Frame_Shift_Del','IGR','Missense_Mutation','Silent','Nonsense_Mutation',
                           'RNA','Splice_Site','Intron','Frame_Shift_Ins','Nonstop_Mutation','In_Frame_Del','ITD','In_Frame_Ins',
                           'Translation_Start_Site',"Multi_Hit", 'Amp', 'Del')
  }else{
    col = colors
  }

  #Default background gray color.
  bg = "#CCCCCC"
  #New version of complexheatmap complains about '', will replace them with random tesx, xxx
  col = c(col, 'xxx' = bg)


  variant.classes = as.character(unique(unlist(as.list(apply(mat, 2, unique)))))
  variant.classes = unique(unlist(strsplit(x = variant.classes, split = ';', fixed = TRUE)))

  variant.classes = variant.classes[!variant.classes %in% c('xxx')]

  type_col = structure(col[variant.classes], names = names(col[variant.classes]))
  type_col = type_col[!is.na(type_col)]

  type_name = structure(variant.classes, names = variant.classes)

  variant.classes = variant.classes[!variant.classes %in% c('Amp', 'Del')]


  if(!is.null(clinicalFeatures)){
    annotationDat = getClinicalData(x = maf)

    if(length(clinicalFeatures[!clinicalFeatures %in% colnames(annotationDat)]) > 0){
      message('Following columns are missing from annotation slot of MAF. Ignoring them..')
      print(clinicalFeatures[!clinicalFeatures %in% colnames(annotationDat)])
      clinicalFeatures = clinicalFeatures[clinicalFeatures %in% colnames(annotationDat)]
      if(length(clinicalFeatures) == 0){
        message('Make sure at-least one of the values from provided clinicalFeatures are present in clinical slot of MAF. Here are available clinical features from MAF..')
        print(colnames(getClinicalData(maf)))
        stop('Zero annotaions to add!')
      }
    }
    annotation = data.frame(row.names = annotationDat$Tumor_Sample_Barcode ,annotationDat[,clinicalFeatures, with = FALSE])
    annotation = annotation[colnames(mat),, drop = FALSE]

    if(!is.null(annotationColor)){
      bot.anno = HeatmapAnnotation(df = annotation, col = annotationColor, annotation_legend_param = list(title_gp = gpar(fontface = "bold"),
                                                                                                          labels_gp = gpar(fontface = "bold")))
    }else{
      bot.anno = HeatmapAnnotation(annotation, annotation_legend_param = list(title_gp = gpar(fontface = "bold"),
                                                                              labels_gp = gpar(fontface = "bold")))
    }
  }


  #------------------------------------Helper functions to add %, rowbar and colbar----------------------------------------------------
  ##This function adds percent rate
  anno_pct = function(index) {
    n = length(index)
    #pct = apply(mat_origin[rev(index), ], 1, function(x) sum(!grepl("^\\s*$", x))/length(x)) * 100
    pct = apply(numMat[rev(index),], 1, function(x) length(x[x != 0]))/as.numeric(maf@summary[3, summary]) * 100
    pct = paste0(round(pct), "%")
    grid::pushViewport(viewport(xscale = c(0, 1), yscale = c(0.5, n + 0.5)))
    grid::grid.text(pct, x = 1, y = seq_along(index), default.units = "native",
                    just = "right", gp = grid::gpar(fontsize = gfs))
    grid::upViewport()
  }

  ha_pct = ComplexHeatmap::HeatmapAnnotation(pct = anno_pct,
                                             width = grid::grobWidth(grid::textGrob("100%", gp = grid::gpar(fontsize = 10, fontface = "bold"))), which = "row")

  ##Following two funcs add grids
  add_oncoprint = function(type, x, y, width, height) {
    grid::grid.rect(x, y, width - unit(0.5, "mm"),
                    height - grid::unit(1, "mm"), gp = grid::gpar(col = NA, fill = bg))

    for (i in 1:length(variant.classes)) {
      if (any(type %in% variant.classes[i])) {
        grid::grid.rect(x, y, width - unit(0.5, "mm"), height -
                          grid::unit(1, "mm"), gp = grid::gpar(col = NA, fill = type_col[variant.classes[i]]))
      } else if (any(type %in% 'Amp')) {
        grid::grid.rect(x, y, width - unit(0.5, "mm"), height -
                          grid::unit(1, "mm"), gp = grid::gpar(col = NA, fill = bg))
        grid::grid.rect(x, y, width - unit(0.5, "mm"), height -
                          unit(15, 'mm'), gp = grid::gpar(col = NA, fill = type_col['Amp']))
      } else if (any(type %in% 'Del')) {
        grid::grid.rect(x, y, width - unit(0.5, "mm"), height -
                          grid::unit(1, "mm"), gp = grid::gpar(col = NA, fill = bg))
        grid::grid.rect(x, y, width - unit(0.5, "mm"), height - grid::unit(15, "mm")
                        , gp = grid::gpar(col = NA, fill = type_col['Del']))
      }
    }
  }

  add_oncoprint2 = function(type, x, y, width, height) {
    for (i in 1:length(variant.classes)) {
      if (any(type %in% variant.classes[i])) {
        grid::grid.rect(x, y, width - unit(0.5, "mm"), height -
                          grid::unit(1, "mm"), gp = grid::gpar(col = NA, fill = type_col[variant.classes[i]]))
      } else if (any(type %in% 'Amp')) {
        grid::grid.rect(x, y, width - unit(0.5, "mm"), height -
                          unit(15, 'mm'), gp = grid::gpar(col = NA, fill = type_col['Amp']))
      } else if (any(type %in% 'Del')) {

        grid::grid.rect(x, y, width - unit(0.5, "mm"), height - grid::unit(15, "mm")
                        , gp = grid::gpar(col = NA, fill = type_col['Del']))
      }
    }
  }

  #This is the main cel function which is passed to ComplexHeatmap::Hetamap()
  celFun = function(j, i, x, y, width, height, fill) {
    type = mat[i, j]
    if(type != 'xxx'){
      typeList = unlist(strsplit(x = as.character(type), split = ';'))
      if(length(typeList) > 1){
        for(i in 1:length(typeList)){
          add_oncoprint2(typeList[i], x, y, width, height)
        }
      }else{
        for(i in 1:length(typeList)){
          add_oncoprint(typeList[i], x, y, width, height)
        }
      }

    }else{
      add_oncoprint(type, x, y, width, height)
    }
  }

  #----------------------------------------------------------------------------------------


  #------Main Heatmap function

  if(showGenes){
    if(!is.null(clinicalFeatures)){
      ht = ComplexHeatmap::Heatmap(matrix = mat,rect_gp = grid::gpar(type = "none"), cell_fun = celFun,
                                   row_names_gp = grid::gpar(fontsize = gfs, fontface = "bold"), show_column_names = showTumorSampleBarcodes,
                                   show_heatmap_legend = FALSE,
                                   column_title = hmName, column_names_gp = grid::gpar(fontsize = fs), bottom_annotation = bot.anno,
                                   column_title_gp = grid::gpar(fontsize = tfs, fontface = "bold"))

    }else{
      ht = ComplexHeatmap::Heatmap(matrix = mat,rect_gp = grid::gpar(type = "none"), cell_fun = celFun,
                                   row_names_gp = grid::gpar(fontsize = gfs, fontface = "bold"), show_column_names = showTumorSampleBarcodes,
                                   show_heatmap_legend = FALSE,
                                   column_title = hmName, column_names_gp = grid::gpar(fontsize = fs),
                                   column_title_gp = grid::gpar(fontsize = tfs, fontface = "bold"))
    }
  }else{
    if(!is.null(clinicalFeatures)){
      ht = ComplexHeatmap::Heatmap(matrix = mat, rect_gp = grid::gpar(type = "none"), cell_fun = celFun,
                                   show_column_names = showTumorSampleBarcodes, show_heatmap_legend = FALSE,
                                   top_annotation_height = grid::unit(2, "cm"), show_row_names = FALSE,
                                   column_title =  hmName, column_names_gp = grid::gpar(fontsize = fs), bottom_annotation = bot.anno,
                                   column_title_gp = grid::gpar(fontsize = tfs, fontface = "bold"))

    }else{
      ht = ComplexHeatmap::Heatmap(matrix = mat, rect_gp = grid::gpar(type = "none"), cell_fun = celFun,
                                   show_column_names = showTumorSampleBarcodes, show_heatmap_legend = FALSE,
                                   top_annotation_height = grid::unit(2, "cm"), show_row_names = FALSE,
                                   column_title =  hmName, column_names_gp = grid::gpar(fontsize = fs),
                                   column_title_gp = grid::gpar(fontsize = tfs, fontface = "bold"))
    }

  }


  #Adding % to left or right
  if(left){
    ht_list = ha_pct + ht
  }else{
    ht_list = ht + ha_pct
  }

  return(list(hm = ht_list, tn = type_name, tc = type_col))
}


#added to allow the user to specify the order of genes from top to bottom (e.g. to correspond with, for example, the odds ratio based on the output of mafCompare
sortByMutationKeepGeneOrder <- function(numMat, genes){

  geneOrder = genes
  numMat = numMat[as.character(geneOrder[geneOrder %in% rownames(numMat)]),, drop = FALSE]

  numMat[numMat != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
  tnumMat = t(numMat) #transposematrix
  numMat = t(tnumMat[do.call(order, c(as.list(as.data.frame(tnumMat)), decreasing = TRUE)), ]) #sort

  return(numMat)
}

