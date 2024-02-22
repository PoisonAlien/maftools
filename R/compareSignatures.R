#' Compares identified denovo mutational signatures to known COSMIC signatures
#'
#' @description Takes results from \code{\link{extractSignatures}} and compares them known COSMIC signatures. Two COSMIC databases are used for comparisons - "legacy" which includes 30 signaures, and "SBS" - which includes updated/refined 65 signatures
#' @details SBS signature database was obtained from https://www.synapse.org/#!Synapse:syn11738319.7
#'
#' @param nmfRes results from \code{\link{extractSignatures}}
#' @param sig_db can be \code{legacy}, \code{SBS}, \code{SBS_v34}. Default \code{SBS_v34}
#' @param verbose Default TRUE
#' @return list containing cosine smilarities, aetiologies if available, and best match.
#' @seealso \code{\link{trinucleotideMatrix}} \code{\link{extractSignatures}} \code{\link{plotSignatures}}
#' @export
#'
compareSignatures = function(nmfRes, sig_db = "SBS_v34", verbose = TRUE){

  sig_db = match.arg(arg = sig_db, choices = c("legacy", "SBS", "SBS_v34"))

  sigs_db = readRDS(file = system.file('extdata', paste0(sig_db, '_signatures.RDs'), package = 'maftools'))
  sigs = sigs_db$db
  aetiology = sigs_db$aetiology

  w = nmfRes$signatures
  sigs = sigs[rownames(w),]

  #corMat = c()
  coSineMat = c()
  for(i in 1:ncol(w)){
    sig = w[,i]
    coSineMat = rbind(coSineMat, apply(sigs, 2, function(x){
      round(crossprod(sig, x)/sqrt(crossprod(x) * crossprod(sig)), digits = 3) #Estimate cosine similarity against all 30 signatures
    }))
    #corMat = rbind(corMat, apply(sigs, 2, function(x) cor.test(x, sig)$estimate[[1]])) #Calulate correlation coeff.
  }
  #rownames(corMat) = colnames(w)
  rownames(coSineMat) = colnames(w)

  best_matches = lapply(1:nrow(coSineMat), function(i){
    ae = aetiology[names(which(coSineMat[i,] == max(coSineMat[i,]))),]
    max_cor = names(which(coSineMat[i,] == max(coSineMat[i,])))
    list(
      aetiology = ae,
      best_match = paste0(
        "Best match: ",
        max_cor,
        " [cosine-similarity: ",
        max(coSineMat[i, ]),
        "]"
      )
    )
  })
  names(best_matches) = rownames(coSineMat)

  if(verbose){
    message('-Comparing against COSMIC signatures')
    message('------------------------------------')
    for(i in 1:nrow(coSineMat)){
      ae = aetiology[names(which(coSineMat[i,] == max(coSineMat[i,]))),]
      ae = paste0("Aetiology: ", ae, " [cosine-similarity: ", max(coSineMat[i,]), "]")
      message('--Found ',rownames(coSineMat)[i], ' most similar to ', names(which(coSineMat[i,] == max(coSineMat[i,]))),sep='')
      message(paste0('   ', ae))
    }
    message('------------------------------------')
  }


  return(list(cosine_similarities = coSineMat, aetiology_db = aetiology, best_match = best_matches))
}
