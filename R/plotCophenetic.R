#' Draw an elbow plot of cophenetic correlation metric.
#' @details This function draws an elbow plot of cophenetic correlation metric.
#' @param res output from \code{\link{estimateSignatures}}
#' @param bestFit rank to highlight. Default NULL
#' @seealso \code{\link{estimateSignatures}} \code{\link{plotCophenetic}}
#' @export
plotCophenetic <- function(res = NULL, bestFit = NULL) {
  if (is.null(res)) {
    stop("Please provide output from estimateSignatures")
  }
  par(mar = c(3, 4, 2, 1))
  plot(res$nmfSummary$rank, res$nmfSummary$cophenetic, axes = FALSE, pch = 16, col = "#D8B365", cex = 1, xlab = NA, ylab = NA, ylim = range(pretty(res$nmfSummary$cophenetic)))
  axis(side = 1, at = res$nmfSummary$rank, labels = res$nmfSummary$rank, lwd = 1, font = 1, cex.axis = 1)
  lines(x = res$nmfSummary$rank, y = round(res$nmfSummary$cophenetic, digits = 4), lwd = 1)
  points(res$nmfSummary$rank, res$nmfSummary$cophenetic, pch = 16, col = "#D8B365", cex = 1.6)
  axis(side = 2, at = pretty(res$nmfSummary$cophenetic), lwd = 1, font = 1, las = 2, cex = 1, cex.axis = 1)
  if (!is.null(bestFit)) {
    segments(x0 = bestFit, y0 = 0, x1 = bestFit, y1 = res$nmfSummary[rank == bestFit, cophenetic], lwd = 1, lty = 1, col = "maroon")
  }
  title(main = "cophenetic metric", adj = 0, font.main = 4)
}
