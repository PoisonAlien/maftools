#' Plot elbow plot for Cophenetic Metric against range of factorization ranks
#' @details Use this function to plot Cophenetic correlation metric against a range of factorization ranks tried. This plot helps to choose optimal factorization rank.
#' Ideally optimal rank is the one at which Cophenetic metric significantly drops, and after which there is little change in the values (elbow).
#' @param nmfRes Signature results from \code{\link{extractSignatures}}
#' @param bestFit Choice of best fit. Can be "auto", or an integer. Default "auto"
#'
#' @export
plotCopheneticMetric = function(nmfRes, bestFit = "auto"){

  nmf.sum = nmfRes$nmfSummary

  if(is.null(nmf.sum)){
    stop("Missing nmfSummary. This is applicable only when extractSignatures is run with nTry.")
  }

  if(!is.null(bestFit)){
    if(bestFit == "auto"){
      #First point where cophenetic correlation coefficient starts decreasing
      bestFit = nmf.sum[diff < 0, rank][1]
    }else if(is.numeric(bestFit)){
      if(!bestFit %in% nmf.sum$rank){
        stop(paste0("bestFit must be within range ", paste(range(nmf.sum$rank), collapse = "..")))
      }
    }
  }

  par(mar = c(3, 3, 2, 1))
  plot(nmf.sum$rank, nmf.sum$cophenetic, axes = FALSE, pch = 16, col = "#D8B365", cex = 1, xlab = NA, ylab = NA, ylim = range(pretty(nmf.sum$cophenetic)))
  axis(side = 1, at = nmf.sum$rank, labels = nmf.sum$rank, lwd = 1, font = 1, cex.axis = 1)
  lines(x = nmf.sum$rank, y = round(nmf.sum$cophenetic, digits = 4), lwd = 1)
  points(nmf.sum$rank, nmf.sum$cophenetic, pch = 16, col = "#D8B365", cex = 1.6)
  axis(side = 2, at = pretty(nmf.sum$cophenetic), lwd = 1, font = 1, las = 2, cex = 1, cex.axis = 1)
  if(!is.null(bestFit)){
    segments(x0 = bestFit, y0 = 0, x1 = bestFit, y1 = nmf.sum[rank == bestFit, cophenetic], lwd= 1, lty = 1, col = "maroon")
  }
  mtext(text = "Factorization rank", side = 1, line = 2)
  title(main = "cophenetic metric", adj = 0, font.main = 4)
}
