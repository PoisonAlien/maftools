detectCP = function(dat){
  dat.spl = split(dat, f = as.factor(dat$Chromosome)) #Split by chromosome
  dat.spl = dat.spl[sapply(X = dat.spl, FUN = nrow) > 6] #Choose chromsomes with at least 5 mutations

  cps = lapply(X = dat.spl, FUN = function(x){
                x = x[order(Start_Position, decreasing = FALSE)] #order mutations acording to position
                x$distance = c(1, diff(x$Start_Position)) #calculate distance between consecutive mutations; in log10 scale
                x_cpt = changepoint::cpt.mean(data = log10(x$distance),
                                              method = 'PELT', class = TRUE, minseglen = 5) #Run change point analysis
                ydat = c()
                if(length(cpts(x_cpt)) > 0){
                  if(length(changepoint::cpts(object = x_cpt)) %% 2 == 0){ #Id odd no of cp's detect don;t do anything
                    xchanges = x[changepoint::cpts(x_cpt)]
                    for(i in seq(1, nrow(xchanges), by = 2)){
                      y = xchanges[i:(i+1)]
                      ycp = data.table::data.table(Chromosome = y[,Chromosome][1],
                                                   Start_Position = y[,Start_Position][1],
                                                   End_Position = y[,End_Position][2])
                      setkey(x = ycp, Chromosome, Start_Position, End_Position)
                      yol = data.table::foverlaps(x = dat, y = ycp, type = "within", nomatch = 0) #overlap and see how nmut per cp
                      ycp[,nMuts := nrow(yol)]
                      ycp[,Avg_intermutation_dist := round(mean(diff(yol[,i.Start_Position])), digits = 2)]
                      ycp[,Size := End_Position - Start_Position]
                      yoll = yol[,.N,.(con.class, Tumor_Sample_Barcode)]
                      ydat = rbind(ydat,
                                   cbind(ycp,
                            data.table::dcast(data = yoll, Tumor_Sample_Barcode ~ con.class, value.var = 'N')))

                    }
                  } #See if odd number of changes
                }
                ydat
              })
  cps = data.table::rbindlist(l = cps, fill = TRUE)
  write.table(x = cps, file = paste(unique(cps[,Tumor_Sample_Barcode]), 'changePoints.tsv', sep='_'), sep='\t', quote = FALSE, row.names = FALSE)
  cps[,Tumor_Sample_Barcode := NULL]
  message('Change points detected at:')
  print(cps)
  cps
}
