
getOncoPlot = function(maf, genes, removeNonMutated = FALSE, colors = NULL, showGenes = TRUE, left = FALSE, showTumorSampleBarcodes = FALSE, hmName = hmName){

  #-----preprocess matrix
  mat_origin = maf@numericMatrix

  if(ncol(mat_origin) < 2){
    stop('Cannot create oncoplot for single sample. Minimum two sample required ! ')
  }

  if(nrow(mat_origin) <2){
    stop('Minimum two genes required !')
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

    #Sort
    mat[mat != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmat = t(mat)

    if(nrow(mat) == 1){
      g = rownames(mat)
      mat = t(tmat[do.call(order, c(as.list(as.data.frame(tmat)), decreasing = TRUE)), ])
      rownames(mat) = g
    }else{
      mat = t(tmat[do.call(order, c(as.list(as.data.frame(tmat)), decreasing = TRUE)), ])
    }


    mat_origin = maf@oncoMatrix
    char.mat = maf@oncoMatrix
    char.mat = char.mat[rownames(mat),, drop = FALSE]
    char.mat = char.mat[,colnames(mat), drop = FALSE]
    mat = char.mat

    if(length(genes.missing) > 0){
      genes.missing.mat = t(matrix(data = '', ncol = ncol(mat), nrow = length(genes.missing)))
      colnames(genes.missing.mat) = genes.missing
      mat = rbind(mat, t(genes.missing.mat))
      #mat = mat[genes,]

      genes.missing.mat2 = t(matrix(data = '', ncol = ncol(mat_origin), nrow = length(genes.missing)))
      colnames(genes.missing.mat2) = genes.missing
      mat_origin = rbind(mat_origin, t(genes.missing.mat2))
    }

    #final matrix for plotting
    mat = mat[genes , , drop = FALSE]
    mat_origin = mat_origin[genes , , drop = FALSE]


  #New version of complexheatmap complains about '' , replacing them with random strinf xxx
  mat[mat == ''] = 'xxx'

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


  #------------------------------------Helper functions to add %, rowbar and colbar----------------------------------------------------
  ##This function adds percent rate
  anno_pct = function(index) {
    n = length(index)
    pct = apply(mat_origin[rev(index), ], 1, function(x) sum(!grepl("^\\s*$", x))/length(x)) * 100
    pct = paste0(round(pct), "%")
    grid::pushViewport(viewport(xscale = c(0, 1), yscale = c(0.5, n + 0.5)))
    grid::grid.text(pct, x = 1, y = seq_along(index), default.units = "native",
                    just = "right", gp = grid::gpar(fontsize = 10))
    grid::upViewport()
  }

  ha_pct = ComplexHeatmap::HeatmapAnnotation(pct = anno_pct,
                                             width = grid::grobWidth(grid::textGrob("100%", gp = grid::gpar(fontsize = 10))), which = "row")

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
    ht = ComplexHeatmap::Heatmap(matrix = mat,rect_gp = grid::gpar(type = "none"), cell_fun = celFun,
                                 row_names_gp = grid::gpar(fontsize = 10), show_column_names = showTumorSampleBarcodes,
                                 show_heatmap_legend = FALSE,
                                 column_title = hmName)
  }else{
    ht = ComplexHeatmap::Heatmap(matrix = mat, rect_gp = grid::gpar(type = "none"), cell_fun = celFun,
                                 show_column_names = showTumorSampleBarcodes, show_heatmap_legend = FALSE,
                                 top_annotation_height = grid::unit(2, "cm"), show_row_names = FALSE,
                                 column_title =  hmName)
  }


  #Adding % to left or right
  if(left){
    ht_list = ha_pct + ht
  }else{
    ht_list = ht + ha_pct
  }

  return(list(hm = ht_list, tn = type_name, tc = type_col))
}
