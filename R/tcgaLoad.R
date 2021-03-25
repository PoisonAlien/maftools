#' Prints available TCGA datasets
#'
#' @description Prints available TCGA cohorts
#' @examples
#' tcgaAvailable()
#' @export
#' @seealso \code{\link{tcgaLoad}}
tcgaAvailable = function(){
  cohorts = data.table::fread(input = "https://github.com/PoisonAlien/TCGAmutations/blob/master/inst/extdata/cohorts.txt?raw=true",
                              verbose = FALSE,
                              showProgress = FALSE)
  cohorts
}
