#' Loads a TCGA cohort
#'
#' @description Loads the user mentioned TCGA cohorts
#' @param study Study names to load. Use  \code{\link{tcgaAvailable}} to see available options.
#' @param source Source for MAF files. Can be \code{MC3} or \code{Firehose}. Default \code{MC3}. Argument may be abbreviated (M or F)
#' @param repo one of "github" (default) and "gitee".
#' @examples
#' # Loads TCGA LAML cohort (default from MC3 project)
#' tcgaLoad(study = "LAML")
#' # Loads TCGA LAML cohort (from Borad Firehose)
#' tcgaLoad(study = "LAML", source = "Firehose")
#' @return An object of class MAF.
#' @export
#' @seealso \code{\link{tcgaAvailable}}
#' @details The function loads curated and pre-compiled MAF objects from TCGA cohorts. TCGA data are obtained from two sources namely, Broad Firehose repository, and MC3 project.
#' @references Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines Kyle Ellrott, Matthew H. Bailey, Gordon Saksena, et. al. Cell Syst. 2018 Mar 28; 6(3): 271–281.e7.

tcgaLoad <- function(study = NULL, source = c("MC3", "Firehose"), repo = c("github", "gitee")) {
  if (is.null(study)) {
    stop("Please provide a study name.\nUse tcgaAvailable() for available cohorts.")
  }
  source <- match.arg(source)

  repo <- match.arg(repo)
  repo_url <- if (repo == "github") {
    "https://github.com/PoisonAlien/TCGAmutations/raw/master"
  } else {
    "https://gitee.com/ShixiangWang/TCGAmutations/raw/master"
  }

  study <- toupper(x = study)
  cohorts <- data.table::fread(
    input = paste0(repo_url, "/inst/extdata/cohorts.txt"),
    verbose = FALSE,
    showProgress = FALSE
  )
  # cohorts = data.table::fread(file = cohorts)
  cohorts <- cohorts[cohorts$Study_Abbreviation %in% study]

  if (nrow(cohorts) == 0) {
    stop("Could not find requested datasets!\nUse tcga_available() for available cohorts.")
  }
  cohort <- as.character(cohorts$Study_Abbreviation)

  study.dat <- paste0(repo_url, "/inst/extdata/", source, "/", cohort, ".RDs")

  if (source == "MC3") {
    doi <- rep("https://doi.org/10.1016/j.cels.2018.03.002", nrow(cohorts))
  } else {
    doi <- cohorts$Firehose
    doi <- unlist(data.table::tstrsplit(x = doi, spli = "\\[", keep = 2))
    doi <- gsub(pattern = "\\]$", replacement = "", x = doi)
  }

  if (length(study.dat) > 1) {
    mafs <- lapply(seq_along(study.dat), function(i) {
      message("Loading ", cohorts$Study_Abbreviation[i], ". Please cite: ", doi[i], " for reference")
      readRDS(file = url(study.dat[i], "rb"))
    })
    names(mafs) <- cohorts$Study_Abbreviation
  } else {
    message("Loading ", cohorts$Study_Abbreviation, ". Please cite: ", doi, " for reference")
    mafs <- readRDS(file = url(study.dat, "rb"))
  }

  mafs
}
