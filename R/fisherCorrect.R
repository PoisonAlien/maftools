fisherCorrection <- function(fc) {
  # Adding a tiny value to zero/Inf odds ratios (zeros become Inf's when log converted)
  fc <- as.data.frame(fc)
  fc$or <- ifelse(test = is.infinite(fc$or), yes = fc$ci.up, no = fc$or)
  fc$or <- ifelse(test = fc$or == 0, yes = fc$ci.low, no = fc$or)

  fc$ci.low <- ifelse(test = fc$ci.low == 0, yes = fc$or, no = fc$ci.low)
  fc$ci.up <- ifelse(test = is.infinite(fc$ci.up), yes = fc$ci.low, no = fc$ci.up)

  return(data.table::data.table(fc))
}
