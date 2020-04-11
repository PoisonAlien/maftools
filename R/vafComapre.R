#' compare VAF across two cohorts
#' @description Draw boxplot distibution of VAFs across two cohorts
#'
#' @param m1 first \code{\link{MAF}} object. Required.
#' @param m2 second \code{\link{MAF}} object. Required.
#' @param genes specify genes for which plot has to be generated. Default NULL.
#' @param top if \code{genes} is NULL plots top n number of genes. Defaults to 5.
#' @param m1Name optional name for first cohort
#' @param m2Name optional name for second cohort
#' @param vafCol1 manually specify column name for vafs in \code{m1}. Default looks for column 't_vaf'
#' @param vafCol2 manually specify column name for vafs in \code{m2}. Default looks for column 't_vaf'
#' @param cols vector of colors corresponding to \code{m1} and \code{m2} respectivelly.
#' @param sigvals Estimate and add significance stars. Default TRUE.
#' @param nrows Number of rows in the layout. Default NULL - estimated automatically
#' @param ncols Number of genes drawn per row. Default 4
#' @export

vafComapre = function(m1, m2, genes = NULL, top = 5, vafCol1 = NULL, vafCol2 = NULL,
                      m1Name = "M1", m2Name = "M2", cols = c("#2196F3", "#4CAF50"), sigvals = TRUE, nrows = NULL, ncols = NULL){

  if(is.null(genes)){
    genes = unique(c(as.character(getGeneSummary(x = m1)[1:top, Hugo_Symbol]),
              as.character(getGeneSummary(x = m2)[1:top, Hugo_Symbol])))
  }


  dat1 = subsetMaf(maf = m1, genes = genes, includeSyn = FALSE, mafObj = FALSE)
  if(!'t_vaf' %in% colnames(dat1)){
    if(is.null(vafCol1)){
      if(all(c('t_ref_count', 't_alt_count') %in% colnames(dat1))){
        message("t_vaf field is missing, but found t_ref_count & t_alt_count columns. Estimating vaf..")
        dat1[,t_vaf := as.numeric(as.character(t_alt_count))/(as.numeric(as.character(t_ref_count)) + as.numeric(as.character(t_alt_count)))]
      }else{
        print(colnames(dat1))
        stop('t_vaf field is missing. Use vafCol to manually specify vaf column name.')
      }
    }else{
      colnames(dat1)[which(colnames(dat1) == vafCol1)] = 't_vaf'
      dat1[,t_vaf := as.numeric(as.character(t_vaf))]
    }
  }

  if(max(dat1$t_vaf, na.rm = TRUE) > 1){
    dat1$t_vaf = dat1$t_vaf/100
  }

  dat2 = subsetMaf(maf = m2, genes = genes, includeSyn = FALSE, mafObj = FALSE)
  if(!'t_vaf' %in% colnames(dat2)){
    if(is.null(vafCol2)){
      if(all(c('t_ref_count', 't_alt_count') %in% colnames(dat2))){
        message("t_vaf field is missing, but found t_ref_count & t_alt_count columns. Estimating vaf..")
        dat2[,t_vaf := as.numeric(as.character(t_alt_count))/(as.numeric(as.character(t_ref_count)) + as.numeric(as.character(t_alt_count)))]
      }else{
        print(colnames(dat2))
        stop('t_vaf field is missing. Use vafCol to manually specify vaf column name.')
      }
    }else{
      colnames(dat2)[which(colnames(dat2) == vafCol2)] = 't_vaf'
      dat2[,t_vaf := as.numeric(as.character(t_vaf))]
    }
  }

  if(max(dat2$t_vaf, na.rm = TRUE) > 1){
    dat2$t_vaf = dat2$t_vaf/100
  }

  vafdat = list(dat1[,.(Hugo_Symbol, t_vaf)], dat2[,.(Hugo_Symbol, t_vaf)])
  names(vafdat) = c(m1Name, m2Name)
  vafdat = data.table::rbindlist(l = vafdat, use.names = TRUE, fill = TRUE, idcol = "Cohort")
  vafdat[,Cohort := factor(x = Cohort, levels = c(m1Name, m2Name), ordered = TRUE)]

  genes = vafdat[,.N,Hugo_Symbol][,Hugo_Symbol]
  vafdat = split(vafdat, as.factor(as.character(vafdat$Hugo_Symbol)))

  cols = c("#2196F3", "#4CAF50")

  if(length(genes) <= 4){
    optrows = 1
    optcols = length(genes)
    matdat = 1:length(genes)
  }else{
    optrows = ceiling(length(genes)/4)
    optcols = 4
    matdat = 1:(ceiling(length(genes)/4)*4)
  }

  if(is.null(nrows)){
    nrows = optrows
  }

  if(is.null(ncols)){
    ncols = optcols
  }

  lo = layout(mat = matrix(data = matdat, nrow = nrows, ncol = ncols, byrow = TRUE))

  par(mar = c(2, 2.5, 2, 1))
  for(i in seq_len(length(vafdat))){
    x = vafdat[[i]]
    b = boxplot(t_vaf ~ Cohort, data = x, xaxt="n", boxwex=0.5, outline=FALSE, lty=1, lwd = 1, outwex=0,
            staplewex=0, ylim = c(0, 1), axes = FALSE, border = cols, horizontal = FALSE, ylab = NA, col = NA, xlim = c(0, 2.5))
    stripchart(t_vaf ~ Cohort, vertical = TRUE, data = x,
               method = "jitter", add = TRUE, pch = 16, col = cols, cex = 0.5)
    if(i == 1){
      axis(side = 2, at = seq(0, 1, 0.2), las = 1, font = 1, lwd = 1, line = 0)
      axis(side = 1, labels = c(m1Name, m2Name), at = c(1, 2))
      #axis(side = 2, at = 1:length(b$names), labels = b$names, tick = FALSE, las = 2, font = 3, line = -1, cex.axis = gene_fs)
    }
    title(main = names(vafdat)[i], adj = 0, font.main = 3)
    abline(h = seq(0, 1, 0.2), v = 1:length(b$names), col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.25), lty = 2)
    if(sigvals){
      sig = try(expr = t.test(t_vaf ~ Cohort, data = x), silent = TRUE)
      if(is(sig, class2 = "try-error")){
        sig = NA
      }
      if(!is.na(sig)){
        sig = sig$p.value
        if(sig < 0.01){
          #segments(x0 = 1, y0 = 1, x1 = 2, y1 = 1, col = "black")
          text(x = 1.5, y = 1, labels = "***")
        }else{
          #segments(x0 = 1, y0 = 1, x1 = 2, y1 = 1, col = "black")
          text(x = 1.5, y = 1, labels = "ns")
        }
      }
    }
  }
}
