#' Prints available TCGA datasets
#'
#' @description Prints available TCGA cohorts
#' @param repo one of "github" (default) and "gitee" (for users in China).
#' @examples
#' tcgaAvailable()
#' @export
#' @seealso \code{\link{tcgaLoad}}
tcgaAvailable = function(repo = c("github", "gitee")){
  repo <- match.arg(repo)
  repo_url <- if (repo == "github") {
    "https://github.com/PoisonAlien/TCGAmutations/raw/master"
  } else {
    "https://gitee.com/ShixiangWang/TCGAmutations/raw/master"
  }
  cohorts = data.table::fread(
    input = paste0(repo_url, "/inst/extdata/cohorts.txt"),
    verbose = FALSE,
    showProgress = FALSE)
  cohorts
}
