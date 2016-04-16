refineClusters = function(clusters){
  totclusts = unique(clusters[,cluster])
  refinedClusters = c()

  for(i in 1:length(totclusts)){
    clust.dat = clusters[cluster == totclusts[i]]
    if(nrow(clust.dat) > 1){
      clust.dat = clusters[cluster == totclusts[i]]

      #Based on boxplot stats
      outl = boxplot.stats(x = clust.dat[,t_vaf])$out
      clust.dat$cluster = ifelse(test = clust.dat$t_vaf %in% outl, yes = 'outlier', no = clust.dat$cluster)

      #Based on zscore
      #clust.dat$z = scale(clust.dat[,t_vaf], center = TRUE, scale = TRUE)
      #clust.dat$cluster = ifelse(test = abs(clust.dat$z) >= 2.5 , no = clust.dat$cluster, yes = 'outlier')
      #clust.dat$cluster = ifelse(test = abs(clust.dat$z) >= 2.5 , no = clust.dat$cluster, yes = 'outlier')
      #clust.dat[,z := NULL]
    }else{
      clust.dat$cluster = 'outlier'
    }
    refinedClusters = rbind(refinedClusters, clust.dat)
  }

  return(refinedClusters)
}
