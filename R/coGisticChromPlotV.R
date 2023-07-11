#' Read the gscore and broad data from a GISTIC2.0 output folder
#'
#' @param gistic_res_dir  the path to the GISTIC2.0 output folder
#'
#' @return LoadedGisticObj object which contains a GISTIC object and a broad data table
#' @export
#'
#' @examples
#' gistic_res_folder = system.file("extdata",package = "maftools")
#' laml.gistic = yload_gistic(gistic_res_folder)
#'
yload_gistic = function(gistic_res_dir){
  res = list()

  lf = list.files(gistic_res_dir)
  lf = lf[stringr::str_starts(lf,'all_lesions.conf_')]
  conf = stringr::str_sub(lf,18,-5)
  res$conf = as.numeric(conf)

  all.lesions <- file.path(gistic_res_dir, paste0("all_lesions.conf_",conf,".txt"))
  amp.genes <- file.path(gistic_res_dir, paste0("amp_genes.conf_",conf,".txt"))
  del.genes <- file.path(gistic_res_dir, paste0("del_genes.conf_",conf,".txt"))
  scores.gis <- file.path(gistic_res_dir, "scores.gistic")
  scores.broad <- file.path(gistic_res_dir, "broad_significance_results.txt")

  stopifnot(file.exists(all.lesions) & file.exists(amp.genes) & file.exists(del.genes) & file.exists(scores.gis))

  res$gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)
  if (file.exists(scores.broad)){
    res$broad = utils::read.table(scores.broad, sep = "\t", header = T, quote = "", stringsAsFactors = F,
                                  comment.char = "", na.strings = "", )
  }else{
    res$broad = NULL
  }

  class(res) = c('LoadedGisticObj',class(res))

  res
}


#' Co-plot version of gisticChromPlot()
#'
#' @description Use two GISTIC object or/and two MAF obejct to view a vertical arranged version of
#'   Gistic Chromosome plot results on the Amp or Del G-scores.
#'
#' @param gistic1 data will be plotted on the left side, yload_gistic() or readGistic() returned
#'   object
#' @param gistic2 data will be plotted on the right side, yload_gistic() or readGistic() returned
#'   object
#' @param g1Name the title of the left side
#' @param g2Name the title of the right side
#' @param type default 'Amp', c('Amp',"Del"), choose one to plot
#' @param markBands default TRUE, integer of length 1 or 2 or TRUE, mark
#' @param labelGenes if you want to label some genes you are interested along the chromosome, set
#'   TRUE
#' @param y_lims default NULL, auto determine the G-score ranges (x axis lims) of the plot. You can
#'   set it to c(-6,6), the left plot values must be negative and right positve, the abs values of
#'   the y_lims will be annotated as the x-axis tick labels. If you set y_lims to not NULL values,
#'   then the `symmetric` will be auto set to FALSE
#' @param maf1,maf2 if labelGenes==TRUE, you need to provide read.maf() object, the genes mutation
#'   info collected from the maf1 is shown on the left side, while maf2 on the right side. the genes
#'   selected are controled by the mutGenes or mutGenes1 or mutGenes1 parameter, see following.
#' @param mutGenes,mutGenes1,mutGenes2 default NULL, could be NULL, number, or character vector of
#'   gene symbols which match the corresponding  MAF object's Hugo_Symbol column values. mutGenes
#'   controls both sides of the annotation, mutGenes1 controls only left side and corresponding data
#'   is extracted from to maf1, and mutGenes2 controls only right side annotation and corresponding
#'   to maf2. If `NULL`, extract the top 50 mutated genes from maf1 and maf2 seperatedly then
#'   annotate them on the left side (maf1 genes) and right side (maf2 genes). if integer, say N,
#'   only top N genes will be extracted seperately from maf1 and maf2. These two condition leads to
#'   different genes annotated on both sides. If character vector, then the genes have mutated in
#'   maf1 and maf2 will be annotated on both side of the figure which mean the two sides have the
#'   same list of genes. if mutGenes is not NULL and both mutGenes1 and mutGenes1 are NULL, then the
#'   auto set mutGenes1 = mutGenes2 = mutGenes.
#' @param fdrCutOff default 0.05,only items with FDR < fdrCutOff will be colored as Amp or Del (
#'   colored 'Red' or 'Blue'), others will be seen as non-significant events (colored gray)
#' @param symmetric default TRUE, If False, when the gistic1 and gistic2 have different max values
#'   of G-scores, the Chrom (0 point of x axis) will not be in the center of the whole plot, if you
#'   set symmetric==TRUE, then the one with smaller max(G-score) will be stretched larger to make
#'   the 0 of the x axis in the middle which eventually make the plot more symmetric.
#' @param color NULL or a named vector. the color of the G-score lines, default NULL which will set
#'   the color c(Amp = "red", Del = "blue", neutral = 'gray70')
#' @param ref.build default "hg19", c('hg18','hg19','hg38') supported at current.(same as maftools
#'   other functions)
#' @param cytobandOffset default 'auto', the width of the chromosome rects(Y axis at 0 point of X
#'   axis). by default will be 1.5% of the width of x axis.
#' @param txtSize the zoom value of most of the texts
#' @param cytobandTxtSize textsize of the cytoband annotation
#' @param mutGenesTxtSize textsize of the mutGenes annotation
#' @param rugTickSize the rug line width of the cytoband annotation
#'
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' all.lesions <- system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
#' amp.genes <- system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
#' del.genes <- system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
#' scores.gis <- system.file("extdata", "scores.gistic", package = "maftools")
#' laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = #' del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)
#'
#' # or you can use new helper function to quickly read gistic results into new object:
#' gistic_res_folder = system.file("extdata",package = "maftools")
#' laml.gistic2 = yload_gistic(gistic_res_folder)
#' # NOTE:
#' # to quickly show the plots, I have to use the data inside the maftools, so the left and right part of the plot is mirrored, but
#' # in practise you will send different GISTIC and MAF objects to `gistic1`,`gistc2`, `maf1`,`maf2`, the function will also work well.
#'
#' laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
#' laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')
#' laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
#' laml2 = laml
#'
#' # --- plot ---
#' gisticChromPlot2v(gistic1 = laml.gistic, gistic2 = laml.gistic2, type='Del',
#'                    symmetric = TRUE, g1Name = 'TCGA1',
#'                    g2Name = 'TCGA2', maf1 = laml, maf2 = laml2, mutGenes = 30)
#' }
coGisticChromPlotV = function(gistic1 = NULL,
                              gistic2 = NULL,
                              g1Name = "",
                              g2Name = "",

                              type = 'Amp',
                              markBands = TRUE,
                              labelGenes = TRUE,
                              y_lims = NULL,

                              maf1 = NULL,
                              maf2 = NULL,

                              mutGenes = NULL,
                              mutGenes1 = NULL,
                              mutGenes2 = NULL,

                              fdrCutOff = 0.05,
                              symmetric = TRUE,

                              color = NULL,
                              ref.build = "hg19",
                              cytobandOffset = 'auto',
                              txtSize = 0.8,
                              cytobandTxtSize = 1,
                              mutGenesTxtSize = 0.6,
                              rugTickSize = 0.1) {
  fdrCutOff.log10 = -log10(fdrCutOff)

  if ('LoadedGisticObj' %in% class(gistic1)) {
    gistic1 = gistic1$gistic
  }
  stopifnot('GISTIC' %in% class(gistic1))

  if ('LoadedGisticObj' %in% class(gistic2)) {
    gistic2 = gistic2$gistic
  }
  stopifnot('GISTIC' %in% class(gistic2))

  #
  g1 = getCytobandSummary(gistic1)

  g1[, `:=`(Chromosome, sapply(strsplit(
    x = g1$Wide_Peak_Limits,
    split = ":"
  ), "[", 1))]
  g1[, `:=`(loc, sapply(strsplit(
    x = g1$Wide_Peak_Limits, split = ":"
  ),
  "[", 2))]
  g1[, `:=`(Start_Position, sapply(strsplit(x = g1$loc, split = "-"),
                                   "[", 1))]
  g1[, `:=`(End_Position, sapply(strsplit(x = g1$loc, split = "-"),
                                 "[", 2))]

  g1.lin = transformSegments(segmentedData = g1[, .(
    Chromosome,
    Start_Position,
    End_Position,
    qvalues,
    Cytoband,
    Variant_Classification
  )])

  gis1.scores = transformSegments(segmentedData = gistic1@gis.scores,
                                  build = ref.build)

  #
  g2 = getCytobandSummary(gistic2)

  g2[, `:=`(Chromosome, sapply(strsplit(
    x = g2$Wide_Peak_Limits,
    split = ":"
  ), "[", 1))]
  g2[, `:=`(loc, sapply(strsplit(
    x = g2$Wide_Peak_Limits, split = ":"
  ),
  "[", 2))]
  g2[, `:=`(Start_Position, sapply(strsplit(x = g2$loc, split = "-"),
                                   "[", 1))]
  g2[, `:=`(End_Position, sapply(strsplit(x = g2$loc, split = "-"),
                                 "[", 2))]

  g2.lin = transformSegments(segmentedData = g2[, .(
    Chromosome,
    Start_Position,
    End_Position,
    qvalues,
    Cytoband,
    Variant_Classification
  )])

  gis2.scores = transformSegments(segmentedData = gistic2@gis.scores,
                                  build = ref.build)



  # g.lin是用来注释bar的
  g1.lin$loc = 'left'
  g2.lin$loc = 'right'
  g.lin = rbind(g1.lin , g2.lin)

  # gis.scores是用来画bar图的
  gis1.scores$loc = 'left'
  gis2.scores$loc = 'right'
  gis.scores = rbind(gis1.scores, gis2.scores)

  if (is.null(color)) {
    COLORS = c(Amp = "red",
               Del = "blue",
               neutral = 'gray70')
  } else{
    if (is.null(names(color))) {
      if (length(color) == 2) {
        COLORS = c(color, 'gray70')
      } else if (length(color) == 3) {
        COLORS = c(color[1], color[3], color[2])
      }
      names(COLORS) = c('Amp', 'Del', 'gray70')
    } else{
      stopifnot(length(intersect(
        c('Amp', 'Del', 'gray70'), names(color)
      )) == 3)
      COLORS = color
    }
  }

  gis.scores$VC = gis.scores$Variant_Classification
  gis.scores$Variant_Classification = ifelse(
    test = as.numeric(gis.scores$fdr) >
      fdrCutOff.log10,
    yes = gis.scores$Variant_Classification,
    no = "neutral"
  )

  gis.scores$Variant_Classification = factor(gis.scores$Variant_Classification,
                                             levels = c("neutral", "Amp", "Del"))
  gis.scores$value = ifelse(
    test = gis.scores$loc ==
      "left",
    yes = -gis.scores$G_Score,
    no = gis.scores$G_Score
  )

  if (ref.build == "hg19") {
    chr.lens = c(
      249250621,
      243199373,
      198022430,
      191154276,
      180915260,
      171115067,
      159138663,
      146364022,
      141213431,
      135534747,
      135006516,
      133851895,
      115169878,
      107349540,
      102531392,
      90354753,
      81195210,
      78077248,
      59128983,
      63025520,
      48129895,
      51304566,
      155270560,
      59373566
    )
  } else if (ref.build == "hg18") {
    chr.lens = c(
      247249719,
      242951149,
      199501827,
      191273063,
      180857866,
      170899992,
      158821424,
      146274826,
      140273252,
      135374737,
      134452384,
      132349534,
      114142980,
      106368585,
      100338915,
      88827254,
      78774742,
      76117153,
      63811651,
      62435964,
      46944323,
      49691432,
      154913754,
      57772954
    )
  } else if (ref.build == "hg38") {
    chr.lens = c(
      248956422,
      242193529,
      198295559,
      190214555,
      181538259,
      170805979,
      159345973,
      145138636,
      138394717,
      133797422,
      135086622,
      133275309,
      114364328,
      107043718,
      101991189,
      90338345,
      83257441,
      80373285,
      58617616,
      64444167,
      46709983,
      50818468,
      156040895,
      57227415
    )
  } else {
    stop("ref.build can only be hg18, hg19 or hg38")
  }

  chr.lens.cumsum = cumsum(chr.lens)
  nchrs = length(unique(gis.scores$Chromosome))
  chr.labels = c(1:22, "X", "Y")
  chr.tbl = data.table::data.table(chr = chr.labels,
                                   start = c(1,
                                             chr.lens.cumsum[1:length(chr.lens.cumsum) - 1]),
                                   end = chr.lens.cumsum)
  chr.tbl$color = rep(c("black", "white"), length = nrow(chr.tbl))

  # 找寻xy轴的范围
  xlims = c(0, chr.lens.cumsum[length(chr.lens.cumsum)])
  if (is.null(y_lims)) {
    y_lims = pretty(gis.scores[gis.scores$VC == type, ][, value], na.rm = TRUE)
  } else {
    # y_lims_basic = pretty(gis.scores[gis.scores$VC==type,][, value], na.rm = TRUE)
    y_lims = pretty(y_lims, na.rm = TRUE)
    # 解决conflicts de custom_lims vs symmetric
    message('As you set y_lims not NULL, forced symmetric = FALSE')
    symmetric = FALSE
  }
  ylims = range(y_lims)
  # 在左右子图lim不同时，是否强制左右对称（放大小的值和label）
  if (symmetric == TRUE) {
    ratio.l = abs(ylims[2] / ylims[1])
    ratio.r = abs(ylims[1] / ylims[2])
    if (ratio.l >= ratio.r) {
      ylims[1] = ylims[1] * ratio.l
      ratio.r = 1
      y_lims = pretty(ylims, na.rm = TRUE)
      y_lims_labels = y_lims
      y_lims_labels[y_lims_labels < 0]  = signif(y_lims_labels[y_lims_labels <
                                                                 0] / ratio.l, 2)
      y_lims_labels = abs(y_lims_labels)
    } else if (ratio.l < ratio.r) {
      ylims[2] = ylims[2] * ratio.r
      ratio.l = 1
      y_lims = pretty(ylims, na.rm = TRUE)
      y_lims_labels = y_lims
      y_lims_labels[y_lims_labels > 0] = signif(y_lims_labels[y_lims_labels > 0] / ratio.r, 2)
      y_lims_labels = abs(y_lims_labels)
    }
  } else{
    ratio.l = 1
    ratio.r = 1
    y_lims_labels = abs(y_lims)
  }

  if (cytobandOffset == 'auto') {
    cytobandOffset = signif((ylims[2] - ylims[1]) * 0.015, 2)
  }


  gis.scores$ystart = ifelse(
    test = gis.scores$loc ==
      "left",
    yes = -cytobandOffset,
    no = cytobandOffset
  )
  # gis.scores$Variant_Classification = factor(x = as.character(gis.scores$Variant_Classification),
  #     levels = c("neutral", "Amp", "Del"))


  gis.scores.splt = split(gis.scores, as.factor(gis.scores$Variant_Classification))

  # make.custom(10,12)
  ## --- plot 开始绘图 ---

  # 设置子fig的layout
  layout.matrix <- matrix(c(2, 1, 3), nrow = 1, ncol = 3)
  layout(mat = layout.matrix, widths = c(2, 8, 2))

  ## --- 画主图 ---
  # 起plot画布
  par(mar = c(4, 4, 0.2, 4)) # c(bottom, left, top, right)
  plot(
    NA,
    NA,
    type = 'l',
    xlim = ylims,
    ylim = xlims,
    axes = FALSE,
    xlab = NA,
    ylab = NA
  )
  # 画组名title（中上外框）
  center = which(y_lims == 0)
  text(
    y = xlims[2],
    x = y_lims[center - 1],
    labels = g1Name ,
    adj = 0,
    cex = 1,
    pos = 3,
    font = 3,
    xpd = TRUE
  )
  text(
    y = xlims[2],
    x = y_lims[center + 1],
    labels = g2Name ,
    adj = 1,
    cex = 1,
    pos = 3,
    font = 3,
    xpd = TRUE
  )
  #mtext(text = g1Name, at =  ,side = 3, line = 0, cex = 1)
  #mtext(text = g2Name, at = y_lims[(length(y_lims) + 3) / 2] ,side = 3, line = 0, cex = 1)

  # 画主数据线
  for (n in names(gis.scores.splt)) {
    if (n == 'neutral' || n == type) {
      # df.and = data.frame.amp.neutral.del
      df.and = gis.scores.splt[[n]]
      colour = COLORS[n]

      df = df.and[loc == 'left']
      segments(
        x0 = df$ystart,
        y0 = xlims[2] - df$Start_Position_updated,
        x1 = df$value * ratio.l + df$ystart,
        y1 = xlims[2] - df$End_Position_updated,
        col = colour,
        lwd = 1.5
      )
      # lines(x = df$value + df$ystart, y = xlims[2] - df$Start_Position_updated, col=colour,lwd=1.5)
      df = df.and[loc == 'right']
      segments(
        x0 = df$ystart,
        y0 = xlims[2] - df$Start_Position_updated,
        x1 = df$value * ratio.r + df$ystart,
        y1 = xlims[2] - df$End_Position_updated,
        col = colour,
        lwd = 1.5
      )
      # lines(x = df$value + df$ystart, y = xlims[2] - df$Start_Position_updated, col=colour,lwd=1.5)
    }
  }


  # 画下面的刻度轴
  axis(
    side = 1,
    at = y_lims,
    las = 1,
    pos = 0,
    labels = y_lims_labels
  )
  mtext(
    text = "G-Score",
    side = 1,
    line = 0.5,
    cex = 1.2
  )
  # 画外框
  rect(
    xleft = ylims[1],
    ybottom = xlims[1],
    xright = ylims[2],
    ytop = xlims[2]
  )
  # 中间染色体框，和染色体label名字
  rect(
    xleft = -cytobandOffset,
    xright = cytobandOffset,
    ytop = xlims[2] - chr.tbl$start,
    ybottom = xlims[2] - chr.tbl$end,
    col = chr.tbl$color
  )
  text(
    y = apply(xlims[2] - chr.tbl[, 2:3], 1, mean),
    x = 0,
    labels = chr.tbl$chr,
    cex = cytobandTxtSize,
    col = c("white", "black")
  )
  # chr 分界线（灰虚线）
  for (i in 1:(length(chr.tbl$end) - 1)) {
    lines(
      y = xlims[2] - rep(chr.tbl$end[i], 2),
      x = ylims,
      lty = 2,
      col = grDevices::adjustcolor("gray80", 0.25)
    )
  }

  # lines(y=xlims,x=rep(fdrCutOff - cytobandOffset,2),lty=2,col = 'green')
  # lines(y=xlims,x=rep(-fdrCutOff + cytobandOffset,2),lty=2,col = 'green')

  # convert numeric (int) markBands into corresponding char markBands names
  if (length(markBands) == 1 && markBands == TRUE) {
    markBands = c(5, 5)
  }
  if (!is.null(markBands)) {
    ordered.g.lin = g.lin[order(qvalues)][Variant_Classification == type]
    if (all(length(markBands) == 1 & markBands == "all")) {
      markBandsdf = ordered.g.lin
    } else if (length(markBands) == 1 &&
               is.numeric(markBands) && markBands >= 1) {
      if (markBands > nrow(g.lin)) {
        markBands = nrow(g.lin)
      }
      markBands = as.integer(markBands)
      markBandsdf = ordered.g.lin[1:markBands, ]
    } else if (length(markBands) == 2 && is.numeric(markBands)) {
      markBands = as.integer(markBands)
      m1 = head(ordered.g.lin[loc == 'left'], markBands[1])
      m2 = head(ordered.g.lin[loc == 'right'], markBands[2])
      markBandsdf = rbind(m1, m2)
    }


    if (nrow(markBandsdf) == 0) {
      message("Available cytobands: ")
      print(getCytobandSummary(x = gistic)[qvalues < fdrCutOff])
      stop(paste(
        "Could not find provided cytobands:",
        paste(markBands,
              collapse = ", ")
      ))
    }

    # gis.scores的区间比markBandsdf(g.lin)的区间要小，"Start_Position_updated", "End_Position_updated"
    # 在两表中不是一一对应, 所以要用foverlap找区间对应关系，而用left_join等keyjoin的方法则不行
    gs = gis.scores.splt[[type]]

    data.table::setkey(x = gs, Start_Position_updated,
                       End_Position_updated)
    cps.l =  data.table::foverlaps(
      by.x = c("Start_Position_updated", "End_Position_updated"),
      by.y = c("Start_Position_updated", "End_Position_updated"),
      x = markBandsdf[loc == 'left'][, .(
        Variant_Classification,
        Cytoband,
        Chromosome,
        Start_Position_updated,
        End_Position_updated,
        loc
      )],
      y = gs[loc == 'left'][, .(
        Start_Position_updated,
        End_Position_updated,
        loc,
        Variant_Classification,
        value,
        ystart,
        Chromosome
      )]
    )
    cps.r = data.table::foverlaps(
      by.x = c("Start_Position_updated", "End_Position_updated"),
      by.y = c("Start_Position_updated", "End_Position_updated"),
      x = markBandsdf[loc == 'right'][, .(
        Variant_Classification,
        Cytoband,
        Chromosome,
        Start_Position_updated,
        End_Position_updated,
        loc
      )],
      y = gs[loc == 'right'][, .(
        Start_Position_updated,
        End_Position_updated,
        loc,
        Variant_Classification,
        value,
        ystart,
        Chromosome
      )]
    )
    #%>% dplyr::arrange(loc,Cytoband,desc(abs(value))) %>% dplyr::distinct(loc,Cytoband,.keep_all = TRUE)
    cyto_peaks_scores = rbind(cps.l, cps.r)
    cyto_peaks_scores = cyto_peaks_scores[order(loc, Cytoband, -abs(value))]
    cyto_peaks_scores = unique(cyto_peaks_scores, by = c("loc", "Cytoband"))
    # drop rows with NA values
    cyto_peaks_scores = cyto_peaks_scores[complete.cases(cyto_peaks_scores)]
    # set label y pos to range(y_lims)
    cyto_peaks_scores$ylims =  signif(ylims[sign(cyto_peaks_scores$ystart) / 2 + 1.5] * 0.9, digits = 2)
    # pos (1=bottom, 2=left, 3=top, 4=right).
    cyto_peaks_scores$pos = sign(cyto_peaks_scores$ystart) + 3

    for (i in 1:nrow(cyto_peaks_scores)) {
      pos = cyto_peaks_scores[i, Start_Position_updated]
      mtext(
        at = xlims[2] - pos,
        # y = cyto_peaks_scores$ylims,#cyto_peaks_scores$amp + cyto_peaks_scores$ystart,
        text = cyto_peaks_scores[i, Cytoband],
        side = cyto_peaks_scores[i, pos],
        font = 3,
        cex = txtSize,
        las = 1
      )
      #     lines(x=p,y=c(1,0.8,0.55),lwd = 0.5, col= mut_dat[i,colour])
      rug(
        x = xlims[2] - pos ,
        side = cyto_peaks_scores[i, pos],
        col = COLORS[type],
        ticksize = rugTickSize * 0.1
      )
    }
  }


  if (!is.null(mutGenes) &&
      is.null(mutGenes1) && is.null(mutGenes2)) {
    if (is.numeric(mutGenes)) {
      mutGenes = as.integer(mutGenes)
    }
    mutGenes1 = mutGenes
    mutGenes2 = mutGenes
    # print(stringr::str_flatten(c(
    #   'both 1,2 use mutGenes values', mutGenes
    # )))
  }

  # 画左Label
  if (!is.null(maf1) && labelGenes) {
    if (is.null(mutGenes1)) {
      mutGenes.app = getGeneSummary(maf1)
      mutGenes.app = head(mutGenes.app, 50)
      mutGenes.app = mutGenes.app[, Hugo_Symbol]
      print('maf1 set, mutGenes.app is null, by default label top 50 genes')
    } else if (is.numeric(mutGenes1)) {
      mutGenes.app = getGeneSummary(maf1)
      mutGenes.app = head(mutGenes.app, mutGenes1)
      mutGenes.app = mutGenes.app[, Hugo_Symbol]
    } else if (is.character(mutGenes1)) {
      mutGenes.app = mutGenes1
    } else{
      stop('when maf1 is set, mutGenes.app must be null or integer')
    }

    mut_dat = transformSegments(segmentedData = maf1@data[,
                                                          .(Chromosome, Start_Position, End_Position, Hugo_Symbol)],
                                build = ref.build)

    # 同gis.scores的区间比markBandsdf, mut_dat(g.lin)的区间比gis1.scores的区间要小，"Start_Position_updated", "End_Position_updated"
    # 在两表中不是一一对应, 所以要用foverlap找区间对应关系，而用left_join等keyjoin的方法则不行
    mut_dat = mut_dat[Hugo_Symbol %chin% mutGenes.app]
    data.table::setkey(x = mut_dat, Start_Position_updated, End_Position_updated)
    gs = gis.scores[VC == type & loc == 'left']
    data.table::setkey(x = gs, Start_Position_updated, End_Position_updated)
    mut_dat = data.table::foverlaps(y = gs, x = mut_dat, mult = "all")

    if (nrow(mut_dat[duplicated(Hugo_Symbol)]) > 0) {
      warning(
        "Multiple CNV region overlaps found for follwing genes. Using the most significant entry for highlighting.",
        immediate. = TRUE
      )
      mut_dat = mut_dat[order(G_Score, decreasing = TRUE)][!duplicated(Hugo_Symbol)]
      mut_dat = mut_dat[complete.cases(mut_dat)]
      #         dups = mut_dat[duplicated(Hugo_Symbol)][, .N, Hugo_Symbol][, Hugo_Symbol]
      #         err = mut_dat[order(Hugo_Symbol, -G_Score)][Hugo_Symbol %in%
      #             dups, .(Hugo_Symbol, Chromosome, Start_Position,
      #             End_Position, Variant_Classification, fdr, G_Score)]
    } else{
      err = NULL
    }

    mut_dat = mut_dat[!is.na(Variant_Classification)]
    mut_dat = mut_dat[order(Chromosome, Start_Position_updated)]
    mut_dat$anno_Position =  round(seq(
      from = xlims[1],
      to = xlims[2],
      length.out = nrow(mut_dat)
    ), 0)

    CLS = COLORS
    CLS['neutral'] = 'black'
    mut_dat$color = CLS[as.character(mut_dat$Variant_Classification)]

    # 绘图
    par(mar = c(4, 0.2, 0.2, 0.2)) # c(bottom, left, top, right)
    # 建坐标轴
    plot(
      NA,
      NA,
      ylim = c(0, chr.lens.cumsum[length(chr.lens.cumsum)]),
      xlim = c(0, 1),
      axes = FALSE,
      xlab = NA,
      ylab = NA
    )
    # 画外框
    # rect(xleft = ylims[1],ybottom = xlims[1],xright = ylims[2],ytop = xlims[2])
    # pos (1=bottom, 2=left, 3=top, 4=right).
    if (nrow(mut_dat) == 0) {
      warning("Could not find mutations")
    } else {
      # 加字
      text(
        y = xlims[2] - mut_dat$anno_Position,
        x = 0.5,
        labels = mut_dat$Hugo_Symbol,
        adj = 2,
        cex = mutGenesTxtSize,
        pos = 2,
        # srt=90,
        font = 3,
        xpd = TRUE,
        col = mut_dat$color
      )
      # 加L型线,x,y按 右→中→左 的顺序绘制
      for (i in 1:nrow(mut_dat)) {
        p = mut_dat[i, c(Start_Position_updated,
                         Start_Position_updated,
                         anno_Position)]
        lines(
          y = xlims[2] - p,
          x = c(1, 0.8, 0.55),
          lwd = 0.5,
          col = mut_dat[i, color]
        )
      }
      # 上下界
      #lines(y=rep(xlims[2],3), x = c(0,0.2,0.45), lwd=2,col='blue')
      #lines(y=rep(0,3), x = c(0,0.2,0.45), lwd=2,col='blue')
    }
  }

  # 画右Label
  if (!is.null(maf2) && labelGenes) {
    if (is.null(mutGenes2)) {
      mutGenes.app = getGeneSummary(maf2)
      mutGenes.app = head(mutGenes.app, 50)
      mutGenes.app = mutGenes.app[, Hugo_Symbol]
      print('maf2 set, mutGenes.app is null, by default label top 50 genes')
    } else if (is.numeric(mutGenes2)) {
      mutGenes.app = getGeneSummary(maf2)
      mutGenes.app = head(mutGenes.app, mutGenes2)
      mutGenes.app = mutGenes.app[, Hugo_Symbol]
    } else if (is.character(mutGenes2)) {
      mutGenes.app = mutGenes2
    } else{
      stop('when maf2 is set, mutGenes.app must be null or integer')
    }

    mut_dat = transformSegments(segmentedData = maf2@data[,
                                                          .(Chromosome, Start_Position, End_Position, Hugo_Symbol)],
                                build = ref.build)

    # 同gis.scores的区间比markBandsdf, mut_dat(g.lin)的区间比gis1.scores的区间要小，"Start_Position_updated", "End_Position_updated"
    # 在两表中不是一一对应, 所以要用foverlap找区间对应关系，而用left_join等keyjoin的方法则不行
    mut_dat = mut_dat[Hugo_Symbol %chin% mutGenes.app]
    data.table::setkey(x = mut_dat, Start_Position_updated, End_Position_updated)
    gs = gis.scores[VC == type & loc == 'right']
    data.table::setkey(x = gs, Start_Position_updated, End_Position_updated)
    mut_dat = data.table::foverlaps(y = gs, x = mut_dat, mult = "all")

    if (nrow(mut_dat[duplicated(Hugo_Symbol)]) > 0) {
      warning(
        "Multiple CNV region overlaps found for follwing genes. Using the most significant entry for highlighting.",
        immediate. = TRUE
      )
      mut_dat = mut_dat[order(G_Score, decreasing = TRUE)][!duplicated(Hugo_Symbol)]
      mut_dat = mut_dat[complete.cases(mut_dat)]
      #         dups = mut_dat[duplicated(Hugo_Symbol)][, .N, Hugo_Symbol][, Hugo_Symbol]
      #         err = mut_dat[order(Hugo_Symbol, -G_Score)][Hugo_Symbol %in%
      #             dups, .(Hugo_Symbol, Chromosome, Start_Position,
      #             End_Position, Variant_Classification, fdr, G_Score)]
    } else{
      err = NULL
    }

    mut_dat = mut_dat[!is.na(Variant_Classification)]
    mut_dat = mut_dat[order(Chromosome, Start_Position_updated)]
    mut_dat$anno_Position = round(seq(
      from = xlims[1],
      to = xlims[2],
      length.out = nrow(mut_dat)
    ), 0)

    CLS = COLORS
    CLS['neutral'] = 'black'
    mut_dat$color = CLS[as.character(mut_dat$Variant_Classification)]

    # 绘图
    par(mar = c(4, 0.2, 0.2, 0.2)) # c(bottom, left, top, right)
    plot(
      NA,
      NA,
      ylim = c(0, chr.lens.cumsum[length(chr.lens.cumsum)]),
      xlim = c(0, 1),
      axes = FALSE,
      xlab = NA,
      ylab = NA
    )
    # 画外框
    # rect(xleft = ylims[1], ybottom = xlims[1], xright = ylims[2],ytop = xlims[2])
    # pos (1=bottom, 2=left, 3=top, 4=right).
    if (nrow(mut_dat) == 0) {
      warning("Could not find mutations")
    } else {
      # 加字
      text(
        y = xlims[2] - mut_dat$anno_Position,
        x = 0.5,
        labels = mut_dat$Hugo_Symbol,
        adj = 2,
        cex = mutGenesTxtSize,
        pos = 4,
        # srt=90,
        font = 3,
        xpd = TRUE,
        col = mut_dat$color
      )
      # 加L型线,x,y按 左→中→右 的顺序绘制
      for (i in 1:nrow(mut_dat)) {
        p = mut_dat[i, c(Start_Position_updated,
                         Start_Position_updated,
                         anno_Position)]
        lines(
          y = xlims[2] - p,
          x = c(0, 0.2, 0.45),
          lwd = 0.5,
          col = mut_dat[i, color]
        )
      }
      # 上下界
      #lines(y=rep(xlims[2],3), x = c(0,0.2,0.45), lwd=2,col='blue')
      #lines(y=rep(0,3), x = c(0,0.2,0.45), lwd=2,col='blue')
    }
  }
  ## only for test usage
  # list(
  #     g1Name = g1Name,
  #     g2Name = g2Name,

  #     type = type,
  #     markBands = markBands,
  #     labelGenes = labelGenes,
  #     y_lims = y_lims,

  #     mutGenes = mutGenes,
  #     mutGenes1 = mutGenes1,
  #     mutGenes2 = mutGenes2,

  #     fdrCutOff = fdrCutOff,
  #     symmetric = symmetric,

  #     color = color,
  #     ref.build = ref.build,
  #     cytobandOffset = cytobandOffset,
  #     txtSize = txtSize,
  #     cytobandTxtSize = cytobandTxtSize,
  #     mutGenesTxtSize = mutGenesTxtSize,
  #     rugTickSize = rugTickSize)
}
