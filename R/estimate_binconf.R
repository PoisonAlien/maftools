#X = number of successes
#n = number of trials
#alpha = pvalues threshold
#Details: https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Agresti-Coull_interval

estimate_binconf = function(X, n, alpha = 0.05){

  z = 1-(alpha/2)
  zsqrd = z^2

  nt = n + zsqrd
  #pr = X/n
  pr = (1/nt) * (X + (zsqrd/2))

  pr_ci_up = pr + z*(sqrt((pr/nt) * (1-pr)))
  pr_ci_low = pr - z*(sqrt((pr/nt) * (1-pr)))

  data.table::data.table(ci_up = pr_ci_up, ci_low = pr_ci_low)
}


#From https://raw.githubusercontent.com/harrelfe/Hmisc/master/R/binconf.s
binconf <- function(x, n, alpha = 0.05,
                    method = c("wilson","exact","asymptotic","all"),
                    include.x = FALSE, include.n = FALSE,
                    return.df = FALSE){
  ## ..modifications for printing and the addition of a
  ##   method argument and the asymptotic interval
  ##   and to accept vector arguments were
  ##   made by Brad Biggerstaff on 10 June 1999

  method <- match.arg(method)
  bc <- function(x, n, alpha, method)
  {
    nu1 <- 2 * (n - x + 1)
    nu2 <- 2 * x
    ll <- if(x > 0)
      x/(x + qf(1 - alpha/2, nu1, nu2) * (n - x + 1))
    else
      0

    nu1p <- nu2 + 2
    nu2p <- nu1 - 2
    pp <- if(x < n)
      qf(1 - alpha/2, nu1p, nu2p)
    else
      1

    ul <- ((x + 1) * pp)/(n - x + (x + 1) * pp)
    zcrit <-  - qnorm(alpha/2)
    z2 <- zcrit * zcrit
    p <- x/n
    cl <- (p + z2/2/n + c(-1, 1) * zcrit *
             sqrt((p * (1 - p) + z2/4/n)/n))/(1 + z2/n)

    if(x == 1)
      cl[1] <-  - log(1 - alpha)/n

    if(x == (n - 1))
      cl[2] <- 1 + log(1 - alpha)/n

    asymp.lcl <- x/n - qnorm(1 - alpha/2) *
      sqrt(((x/n) * (1 - x/n))/n)

    asymp.ucl <- x/n + qnorm(1 - alpha/2) * sqrt(((x/n) * (1 - x/n)
    )/n)
    res <- rbind(c(ll, ul), cl, c(asymp.lcl, asymp.ucl))
    res <- cbind(rep(x/n, 3), res)

    ##dimnames(res) <- list(c("Exact", "Wilson", "Asymptotic"), c(
    ## "Point Estimate", "Lower", "Upper"))
    switch(method,
           wilson =     res[2,  ],
           exact =      res[1,  ],
           asymptotic = res[3,  ],
           all =        res,
           res)
  }

  if((length(x) != length(n)) & length(x) == 1)
    x <- rep(x, length(n))
  if((length(x) != length(n)) & length(n) == 1)
    n <- rep(n, length(x))
  if((length(x) > 1 | length(n) > 1) & method == "all") {
    method <- "wilson"
    warning("method=all will not work with vectors...setting method to wilson")
  }
  if(method == "all" & length(x) == 1 & length(n) == 1) {
    mat <- bc(x, n, alpha, method)
    dimnames(mat) <- list(c("Exact", "Wilson", "Asymptotic"),
                          c("PointEst", "Lower", "Upper"))
    if(include.n)
      mat <- cbind(N = n, mat)

    if(include.x)
      mat <- cbind(X = x, mat)

    if(return.df)
      mat <- as.data.frame(mat)

    return(mat)
  }

  mat <- matrix(ncol = 3, nrow = length(x))
  for(i in 1:length(x))
    mat[i,  ] <- bc(x[i], n[i], alpha = alpha, method = method)

  dimnames(mat) <- list(rep("", dim(mat)[1]),
                        c("PointEst", "Lower", "Upper"))
  if(include.n)
    mat <- cbind(N = n, mat)

  if(include.x)
    mat <- cbind(X = x, mat)

  if(return.df)
    mat <- as.data.frame(mat, row.names=NULL)

  mat
}
