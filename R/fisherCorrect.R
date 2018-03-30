fisherCorrection = function(fc){
   #Adding a tiny value to zero/Inf odds ratios (zeros become Inf's when log converted)
  fc = as.data.frame(fc)
  fc$or = ifelse(test = is.infinite(fc$or), yes = fc$ci.up, no = fc$or )
  fc$or = ifelse(test = fc$or == 0, yes = fc$ci.low, no = fc$or)

  fc$ci.up = ifelse(test = fc$ci.up == 0, yes = fc$or, no = fc$ci.up)
  fc$ci.low = ifelse(test = is.infinite(fc$ci.low), yes = fc$ci.up, no = fc$ci.low)

  return(data.table::data.table(fc))

}
