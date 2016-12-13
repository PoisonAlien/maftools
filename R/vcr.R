#' Samll internal function to make complex events.
#' @description Samll internal function to make complex events. Ignore this.
#' @param xstr charcter to split
#' @param gis Is input from gistic. Logical.
#' @return split string
#' @export

vcr = function(xstr, gis = FALSE) {
  x = as.character(xstr)
  x = strsplit(x = x, split = ';', fixed = TRUE)[[1]]
  x = unique(x)
  xad = x[x %in% c('Amp', 'Del')]
  xvc = x[!x %in% c('Amp', 'Del')]

  if(gis){
    x = ifelse(test = length(xad) > 1, no = xad, yes = 'Complex')
  }else{
    if(length(xvc)>0){
      xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
    }
    if(length(xad) > 0){
      x = paste(xad, xvc, sep = ';')
    }else{
      x = xvc
    }
  }
  return(x)
}
