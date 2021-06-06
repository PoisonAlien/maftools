#' Draw two barplots side by side for cohort comparision.
#' @details Draws two barplots side by side to display difference between two cohorts.
#' @param m1 first \code{\link{MAF}} object
#' @param m2 second \code{\link{MAF}} object
#' @param genes genes to be drawn. Default takes top 5 mutated genes.
#' @param orderBy Order genes by mutation rate in `m1` or `m2`. Default `NULL`, keeps the same order of `genes`
#' @param m1Name optional name for first cohort
#' @param m2Name optional name for second cohort
#' @param colors named vector of colors for each Variant_Classification.
#' @param normalize Default TRUE.
#' @param yLims Default NULL. Auto estimates. Maximum values for `m1` and `m2` respectively
#' @param borderCol Default gray
#' @param titleSize Default 1
#' @param geneSize Default 0.8
#' @param showPct Default TRUE
#' @param pctSize Default 0.7
#' @param axisSize Default 0.8
#' @param showLegend Default TRUE.
#' @param legendTxtSize Default 0.8
#' @param geneMar Default 4
#' @export
#' @examples
#' #' ##Primary and Relapse APL
#' primary.apl <- system.file("extdata", "APL_primary.maf.gz", package = "maftools")
#' relapse.apl <- system.file("extdata", "APL_relapse.maf.gz", package = "maftools")
#' ##Read mafs
#' primary.apl <- read.maf(maf = primary.apl)
#' relapse.apl <- read.maf(maf = relapse.apl)
#' ##Plot
#' coBarplot(m1 = primary.apl, m2 = relapse.apl, m1Name = 'Primary APL', m2Name = 'Relapse APL')
#' dev.off()
#' @return Returns nothing. Just draws plot.
#'
coBarplot = function(m1, m2, genes = NULL, orderBy = NULL,
                     m1Name = NULL, m2Name = NULL, colors = NULL,
                     normalize = TRUE, yLims = NULL, borderCol = "gray",
                     titleSize = 1, geneSize = 0.8,
                     showPct = TRUE, pctSize = 0.7, axisSize = 0.8,
                     showLegend = TRUE, legendTxtSize = 1, geneMar = 4){

  if(is.null(genes)){
    m1.genes = getGeneSummary(m1)[1:5, Hugo_Symbol]
    m2.genes = getGeneSummary(m2)[1:5, Hugo_Symbol]
    genes = rev(unique(c(m1.genes, m2.genes)))
  }

  if(is.null(colors)){
    colors = get_vcColors()
  }

  m1.gs = get_col_df(m = m1, g = genes)
  m1.ss = as.numeric(m1@summary[ID %in% "Samples", summary])
  gs1.load = getGeneSummary(x = m1)[Hugo_Symbol %in% genes,.(Hugo_Symbol, AlteredSamples)]
  gs1.load[,AlteredSamples := round(AlteredSamples/m1.ss, digits = 2) * 100]
  data.table::setDF(x = gs1.load, rownames = gs1.load$Hugo_Symbol)
  gs1.load = gs1.load[, "AlteredSamples", drop = FALSE]
  m1.missing = genes[!genes %in% rownames(gs1.load)]
  if(length(m1.missing) > 0){
    gs1.load = rbind(gs1.load, data.frame(row.names = m1.missing, AlteredSamples = rep(0, length(m1.missing)), stringsAsFactors = FALSE))
  }

  m2.gs = get_col_df(m = m2, g = genes)

  if (ncol(m1.gs) == 0L & ncol(m2.gs) == 0L) {
    stop("Cannot plot genes without mutation!")
  } else if (ncol(m1.gs) == 0L) {
    m1.gs <- m2.gs
    m1.gs[] <- 0
  } else if (ncol(m2.gs) == 0L) {
    m2.gs <- m1.gs
    m2.gs[] <- 0
  }

  m2.ss = as.numeric(m2@summary[ID %in% "Samples", summary])
  gs2.load = getGeneSummary(x = m2)[Hugo_Symbol %in% genes,.(Hugo_Symbol, AlteredSamples)]
  gs2.load[,AlteredSamples := round(AlteredSamples/m2.ss, digits = 2) * 100]
  data.table::setDF(x = gs2.load, rownames = gs2.load$Hugo_Symbol)
  gs2.load = gs2.load[, "AlteredSamples", drop = FALSE]
  m2.missing = genes[!genes %in% rownames(gs2.load)]
  if(length(m2.missing) > 0){
    gs2.load = rbind(gs2.load, data.frame(row.names = m2.missing, AlteredSamples = rep(0, length(m2.missing)), stringsAsFactors = FALSE))
  }

  if(!is.null(orderBy)){
    orderBy = match.arg(arg = orderBy, choices = c("m1", "m2"), several.ok = FALSE)
    if(orderBy == "m1"){
      genes = rev(rownames(gs1.load[order(gs1.load$AlteredSamples, decreasing = TRUE),, drop = FALSE]))
    }else if(orderBy == "m2"){
      genes = rev(rownames(gs2.load[order(gs2.load$AlteredSamples, decreasing = TRUE),, drop = FALSE]))
    }
  }

  m1.gs = t(m1.gs[genes,,drop = FALSE])
  gs1.load = gs1.load[genes, , drop = FALSE]
  m2.gs = t(m2.gs[genes,,drop = FALSE])
  gs2.load = gs2.load[genes, , drop = FALSE]

  if(normalize){
    m1.gs = (m1.gs/m1.ss)*100
    m2.gs = (m2.gs/m2.ss)*100
    if(!is.null(yLims)){
      xat = pretty(c(-yLims[1], 0, yLims[2]))
    }else{
      xat = seq(-100, 100, 20)
    }
  }else{
    if(!is.null(yLims)){
      xat = pretty(c(-yLims[1], 0, yLims[2]))
    }else{
      m1max = max(apply(m1.gs, 2, function(x) cumsum(x)))
      m2max = max(apply(m2.gs, 2, function(x) cumsum(x)))
      m12max = max(m1max, m2max)
      xat = pretty(c(-m12max, 0, m12max))
    }
  }

  lo = matrix(data = c(1, 1, 2, 2), nrow = 2, byrow = TRUE)
  if (showLegend) graphics::layout(mat = lo, heights = c(4, 1.25))
  else graphics::layout(mat = lo, heights = c(20, 1.25))

  par(mar = c(2, geneMar, 1, 2))
  b1 = barplot(
    height = -m1.gs,
    horiz = TRUE,
    las = TRUE,
    col = colors[rownames(m1.gs)],
    xlim = range(xat),
    axes = FALSE, names.arg = rep(NA, ncol(m1.gs)), border = borderCol
  )

  if(showPct){
    text(x = -colSums(m1.gs), y = b1, labels = paste0(gs1.load$AlteredSamples, "%"), adj = 1.1, cex = pctSize, xpd = TRUE)
  }

  barplot(
    height = m2.gs,
    horiz = TRUE,
    las = TRUE,
    col = colors[rownames(m2.gs)],
    add = TRUE,
    axes = FALSE, names.arg = rep(NA, ncol(m2.gs)), border = borderCol
  )
  if(showPct){
    text(x = colSums(m2.gs), y = b1, labels = paste0(gs2.load$AlteredSamples, "%"), adj = -0.1, cex = pctSize, xpd = TRUE)
  }

  if(normalize){
    axis(side = 1, at = xat, labels = paste0(abs(xat), "%"), cex.axis = axisSize)
  }else{
    axis(side = 1, at = xat, labels = abs(xat), cex.axis = axisSize)
  }

  mtext(text = colnames(m1.gs), side = 2, font = 3, las = 2, at = b1, line = 0.75, cex = geneSize)
  if(normalize){
    mtext(text = "Percent of cases", side = 1, line = 2)
  }else{
    mtext(text = "Number of cases", side = 1, line = 2)
  }

  title(main = paste0(m1Name, ' [N=', m1.ss, ']'), cex.main = titleSize, outer = FALSE, font = 2, adj = 0)
  title(main = paste0(m2Name, ' [N=', m2.ss, ']'), cex.main = titleSize, outer = FALSE, font = 2, adj = 1)

  par(mar = c(0, 0.5, 1, 0), xpd = TRUE)
  plot(NULL,ylab='',xlab='', xlim=0:1, ylim=0:1, axes = FALSE)

  vcs = unique(c(rownames(m1.gs), rownames(m2.gs)))
  col = colors[vcs]

  if (showLegend) {
    legend("topleft", legend = names(col), col = col,  bty = "n", border=NA,
           xpd = TRUE, text.font = 1, pch = 15, xjust = 0, yjust = 0,
           cex = legendTxtSize, y.intersp = 1.5, x.intersp = 1,
           pt.cex = 1.2 * legendTxtSize, ncol = ceiling(length(col)/4))
  }

  invisible(list(m1 = m1.gs, m2 = m2.gs))
}

get_col_df = function(m, g){
  if(nrow(getGeneSummary(x = m)[Hugo_Symbol %in% g]) == 0){
    return(data.frame(row.names = g, stringsAsFactors = FALSE))
  }
  ml = apply(createOncoMatrix(
    m = m, g = g, chatty = FALSE, add_missing = TRUE)[['oncoMatrix']],
    1, table, simplify = FALSE)
  ml = lapply(ml, function(x){
    data.frame(x)
  })
  ml = data.table::rbindlist(l = ml, use.names = TRUE, fill = TRUE, idcol = "Gene")
  ml = ml[!Var1 %in% ""]
  #CNV+Mutated = Complex_Event
  ml$Var1 = ifelse(test = ml$Var1 %like% ";", yes = "Complex_Event", no = as.character(ml$Var1))
  ml = ml[,sum(Freq), .(Gene, Var1)]
  colnames(ml) = c("Gene", "Var1", "Freq")
  ml$Gene = factor(x = ml$Gene, levels = g, ordered = TRUE)
  ml$Var1 = as.character(ml$Var1)
  mdf = data.table::dcast(data = ml, Gene ~ Var1, value.var = "Freq", fill = 0, drop = FALSE)
  data.table::setDF(x = mdf, rownames = as.character(mdf$Gene))
  mdf$Gene = NULL
  mdf
}
