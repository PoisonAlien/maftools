detectCP = function(dat, segLen = 5){

  # x_cpt = changepoint::cpt.mean(data = dat$diff, method = 'PELT', class = TRUE, minseglen = 1) #Run change point analysis
  # xchanges = dat[changepoint::cpts(object = x_cpt)]
  # dat.spl = split(xchanges, as.factor(xchanges$Chromosome))
  #
  # if(nrow(xchanges) == 0){
  #   return(NULL)
  # }
  #
  # cps = lapply(X = dat.spl, FUN = function(xcp){
  #   ydat = c()
  #
  #   if(nrow(xcp) > 0){
  #     if(nrow(xcp) %% 2 == 0){ #Id odd no of cp's detect don;t do anything
  #       for(i in seq(1, nrow(xcp), by = 2)){
  #         y = xcp[i:(i+1)]
  #         ycp = data.table::data.table(Chromosome = y[,Chromosome][1],
  #                                      Start_Position = y[,Start_Position][1],
  #                                      End_Position = y[,End_Position][2])
  #         setkey(x = ycp, Chromosome, Start_Position, End_Position)
  #         yol = data.table::foverlaps(x = dat, y = ycp, type = "within", nomatch = 0) #overlap and see how mut per cp
  #         ycp[,nMuts := nrow(yol)]
  #         ycp[,Avg_intermutation_dist := round(mean(diff(yol[,i.Start_Position])), digits = 2)]
  #         ycp[,Size := End_Position - Start_Position]
  #         yoll = yol[,.N,.(con.class, Tumor_Sample_Barcode)]
  #         ydat = rbind(ydat,
  #                      cbind(ycp,
  #                            data.table::dcast(data = yoll, Tumor_Sample_Barcode ~ con.class, value.var = 'N')), fill = TRUE)
  #
  #       }
  #     } #See if odd number of changes
  #   }
  #   ydat
  # })

  # Earlier method: Run per chormosome
  dat.spl = split(dat, f = as.factor(dat$Chromosome)) #Split by chromosome
  dat.spl = dat.spl[sapply(X = dat.spl, FUN = nrow) > 6] #Choose chromsomes with at least 5 mutations
  cps = lapply(X = dat.spl, FUN = function(x){
                x = x[order(as.numeric(as.character(Start_Position)), decreasing = FALSE)] #order mutations acording to position
                x$distance = c(1, diff(as.numeric(as.character(x$Start_Position)))) #calculate distance between consecutive mutations; in log10 scale

                #x_cpt = changepoint::cpt.mean(data = log10(x$distance), penalty = "Asymptotic", pen.value = 0.01,method = 'AMOC', class = TRUE) #Run change point analysis
                x_cpt = try(expr = changepoint::cpt.mean(data = log10(x$distance), method = 'PELT', class = TRUE, minseglen = segLen, ...), silent = TRUE) #Run change point analysis

                ydat = c()
                if(class(x_cpt) != "try-error"){
                  if(length(changepoint::cpts(x_cpt)) > 0){
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
                                           data.table::dcast(data = yoll, Tumor_Sample_Barcode ~ con.class, value.var = 'N')), fill = TRUE)

                      }
                    } #See if odd number of changes
                  }
                }
                ydat
              })

  cps = data.table::rbindlist(l = cps, fill = TRUE)


  if(nrow(cps) > 0){
    #Use Alexandrov's filters:
    #Putative regions of kataegis were identified as those segments containing six or more consecutive mutations
    #with an average intermutation distance of less than or equal to 1,000 bp.
    cps$Filter = ifelse(test = cps$nMuts >= 6 & cps$Avg_intermutation_dist <= 1000, yes = "PASS", no = ".")

    if(nrow(cps) > 0){
      write.table(x = cps, file = paste(unique(cps[,Tumor_Sample_Barcode]), 'changePoints.tsv', sep='_'), sep='\t', quote = FALSE, row.names = FALSE)
      cps[,Tumor_Sample_Barcode := NULL]
      message('Change points detected at:')
      print(cps)
      return(cps)
    }else{
      return(NULL)
    }
  } else{
    return(NULL)
  }
}
