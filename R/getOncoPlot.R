getOncoPlot = function(maf, genes, removeNonMutated = FALSE, colors = NULL, showGenes = TRUE, left = FALSE, showTumorSampleBarcodes = FALSE, hmName = hmName){

  mat_origin = maf@numericMatrix

  if(ncol(mat_origin) < 2){
    stop('Cannot create oncoplot for single sample. Minimum two sample required ! ')
  }

  if(nrow(mat_origin) <2){
    stop('Minimum two genes required !')
  }

  genes.missing = genes[!genes %in% rownames(mat_origin)]
  genes.present = genes[genes %in% rownames(mat_origin)]
  mat = mat_origin[genes.present,]

  #remove nonmutated samples to improve visualization
  if(removeNonMutated){
    tsb = colnames(mat)
    tsb.exclude = colnames(mat[,colSums(mat) == 0])
    tsb.include = tsb[!tsb %in% tsb.exclude]
    mat = mat[,tsb.include]
  }

  #Sort
  mat[mat != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
  tmat = t(mat)
  mat = t(tmat[do.call(order, c(as.list(as.data.frame(tmat)), decreasing = TRUE)), ])

  mat_origin = maf@oncoMatrix
  char.mat = maf@oncoMatrix
  char.mat = char.mat[rownames(mat),]
  char.mat = char.mat[,colnames(mat)]
  #final matrix for plotting
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

  mat = mat[genes,]
  mat_origin = mat_origin[genes,]

  if(is.null(colors)){
    col = c(RColorBrewer::brewer.pal(12,name = "Paired"),RColorBrewer::brewer.pal(11,name = "Spectral")[1:3],'black')
    names(col) = names = c('Nonstop_Mutation','Frame_Shift_Del','Silent','Missense_Mutation','IGR','Nonsense_Mutation',
                           'RNA','Splice_Site','Intron','Frame_Shift_Ins','In_Frame_Dell','In_Frame_Del','ITD','In_Frame_Ins','Translation_Start_Site',"Multi_Hit")
  }else{
    col = colors
  }


  tc = unique(unlist(apply(mat,1,unique)))
  tc = tc[!tc=='']

  type_col = col[tc]
  type_name = names(type_col)
  names(type_name) = type_name

  #from oncoprint source ComplexHeatmap
  add_oncoprint = function(type, x, y, width, height) {

    for(i in 1:length(type_name)){
      if(any(type %in% type_name[i])) {
        grid::grid.rect(x, y, width - grid::unit(0.5, "mm"), height - grid::unit(1, "mm"), gp = grid::gpar(col = NA, fill = type_col[type_name[i]]))
      }
    }
    if(any(type %in% "")) {
      grid::grid.rect(x, y, width - grid::unit(0.5, "mm"), height - grid::unit(1, "mm"), gp = grid::gpar(col = NA, fill = "#CCCCCC"))
    }
  }

  anno_pct = function(index) {
    n = length(index)
    pct = apply(mat_origin[rev(index), ], 1, function(x) sum(!grepl("^\\s*$", x))/length(x)) * 100
    pct = paste0(round(pct), "%")
    grid::pushViewport(viewport(xscale = c(0, 1), yscale = c(0.5, n + 0.5)))
    grid::grid.text(pct, x = 1, y = seq_along(index), default.units = "native",
                    just = "right", gp = grid::gpar(fontsize = 10))
    grid::upViewport()
  }

  ha_pct = ComplexHeatmap::HeatmapAnnotation(pct = anno_pct, width = grid::grobWidth(grid::textGrob("100%", gp = grid::gpar(fontsize = 10))), which = "row")


  if(showGenes){
    ht = ComplexHeatmap::Heatmap(mat, rect_gp = grid::gpar(type = "none"), cell_fun = function(j, i, x, y, width, height, fill) {
      type = mat[i,j]
      add_oncoprint(type, x, y, width, height)
    }, row_names_gp = grid::gpar(fontsize = 10), show_column_names = showTumorSampleBarcodes, show_heatmap_legend = FALSE,
    top_annotation_height = grid::unit(2, "cm"), column_title = hmName)
  }else{
    ht = ComplexHeatmap::Heatmap(mat, rect_gp = grid::gpar(type = "none"), cell_fun = function(j, i, x, y, width, height, fill) {
      type = mat[i,j]
      add_oncoprint(type, x, y, width, height)
    }, show_column_names = showTumorSampleBarcodes, show_heatmap_legend = FALSE,
    top_annotation_height = grid::unit(2, "cm"), show_row_names = FALSE, column_title =  hmName)
  }

  if(left){
    ht_list = ha_pct + ht
  }else{
    ht_list = ht + ha_pct
  }


  #legend = grid::legendGrob(labels = type_name[names(type_col)],  pch = 15, gp = grid::gpar(col = type_col), nrow = 2)

  return(list(hm = ht_list, tn = type_name, tc = type_col))

}
