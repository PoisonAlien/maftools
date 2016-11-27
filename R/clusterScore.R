#--------------------- based on binaomial distribution, estimate threshhold.
get_threshold = function(gene_muts, gene_length){
  th = which(unlist(lapply(X = 2:gene_muts, FUN = function(x) dbinom(x = x, size = gene_muts, prob = 1/gene_length) )) < 0.01)[1]
  return(th+1)
}
#-------------------- end of function.


#--------------------- parse protein positions from protein conversions.
parse_prot = function(dat, AACol, gl, m, calBg = FALSE, nBg){

  if(is.null(AACol)){
    if(! 'AAChange' %in% colnames(dat)){
      message('Available fields:')
      print(colnames(dat))
      stop('AAChange field not found in MAF. Use argument AACol to manually specifiy field name containing protein changes.')
    }
  }else{
    colnames(dat)[which(colnames(dat) == AACol)] = 'AAChange'
  }

  all.prot.dat = dat[,.(Hugo_Symbol, Variant_Classification, AAChange)]
  all.prot.dat = all.prot.dat[Variant_Classification != 'Splice_Site']
  #parse AAchanges to get postion
  prot.spl = strsplit(x = as.character(all.prot.dat$AAChange), split = '.', fixed = TRUE)
  prot.conv = sapply(prot.spl, function(x) x[length(x)])

  all.prot.dat[,conv := prot.conv]
  all.prot.dat = all.prot.dat[!conv == 'NULL']

  #If conversions are in HGVSp_long (default HGVSp) format, we will remove strings Ter followed by anything (e.g; p.Asn1986GlnfsTer13)
  pos = gsub(pattern = 'Ter.*', replacement = '',x = all.prot.dat$conv)

  #Following parsing takes care of most of HGVSp_short and HGVSp_long format
  pos = gsub(pattern = '[[:alpha:]]', replacement = '', x = pos)
  pos = gsub(pattern = '\\*$', replacement = '', x = pos) #Remove * if nonsense mutation ends with *
  pos = gsub(pattern = '^\\*', replacement = '', x = pos) #Remove * if nonsense mutation starts with *
  pos = gsub(pattern = '\\*.*', replacement = '', x = pos) #Remove * followed by position e.g, p.C229Lfs*18

  pos = suppressWarnings( as.numeric(sapply(strsplit(x = pos, split = '_', fixed = TRUE), '[', 1)) )
  all.prot.dat[,pos := pos]

  if(nrow( all.prot.dat[is.na(all.prot.dat$pos),]) > 0){
    #message(paste('Removed', nrow( all.prot.dat[is.na(all.prot.dat$pos),]), 'mutations for which AA position was not available', sep = ' '))
    #print(prot.dat[is.na(prot.dat$pos),])
    all.prot.dat = all.prot.dat[!is.na(all.prot.dat$pos),]
  }

  gene.sum = summarizeMaf(maf = dat, chatty = FALSE)$gene.summary
  #gene.sum = merge.data.frame(x = gene.sum, y = gl, by = 'Hugo_Symbol', all.x = TRUE)
  gene.sum = merge(x = gene.sum, y = gl, by = 'Hugo_Symbol', all.x = TRUE)
  #gene.sum = gene.sum[!is.na(gene.sum$aa.length),]
  gene.sum = gene.sum[!is.na(gene.sum$aa.length)]

  num_mut_colIndex = which(colnames(gene.sum) == 'total')
  aalen_colIndex = which(colnames(gene.sum) == 'aa.length')

  #Get background threshold
  gene.sum$th = apply(gene.sum, 1, function(x) get_threshold(gene_muts = as.numeric(x[num_mut_colIndex]), gene_length = as.numeric(x[aalen_colIndex])))
  #use only genes with atleast 2 (or m ) mutations.
  gene.sum = gene.sum[total >= m]

  if(calBg){
    if(nrow(gene.sum) < nBg){
      #message("Not enough genes to build background. Using predefined values. (Mean = 0.279; SD = 0.13)")
      return(NULL)
    } else{
      syn.res = c()
      pb <- txtProgressBar(min = 0, max = nrow(gene.sum), style = 3) #progress bar

      for(i in 1:nrow(gene.sum)){
        prot.dat = all.prot.dat[Hugo_Symbol == gene.sum[i, "Hugo_Symbol"]]
        syn.res = rbind(syn.res, cluster_prot(prot.dat = prot.dat, gene = gene.sum[i, "Hugo_Symbol"], th = gene.sum[i,"th"], protLen = gene.sum[i,"aa.length"]))
        setTxtProgressBar(pb, i)
      }
      return(syn.res)
    }
  } else{
    nonsyn.res = c()
    pb <- txtProgressBar(min = 0, max = nrow(gene.sum), style = 3) #progress bar

    for(i in 1:nrow(gene.sum)){
      hs = gene.sum[i, Hugo_Symbol]
      #print(hs)
      prot.dat = all.prot.dat[Hugo_Symbol %in% hs]
      nonsyn.res = rbind(nonsyn.res, cluster_prot(prot.dat = prot.dat, gene = hs, th = gene.sum[Hugo_Symbol %in% hs, th], protLen = gene.sum[Hugo_Symbol %in% hs, aa.length]))
      setTxtProgressBar(pb, i)
    }
    return(nonsyn.res)
  }
}

#--------------------- End of Function


# -------------------Clustering function-------------------
cluster_prot = function(prot.dat, gene, th, protLen){

  mergeDist = 5 #hard coded inter event distance.
  #prot.dat = all.prot.dat[Hugo_Symbol == gene]

  #Summarise counts per position
  pos.counts = prot.dat[,.N,pos]
  pos.counts = pos.counts[order(pos)]

  #classify position as meaningful if its greater than background threshhold.
  pos.counts$cluster = ifelse(test = pos.counts$N >= th, yes = 'meaningful', no = 'nonMeaningful')

  #Just choose meaningful positions
  clust.tbl = pos.counts[cluster %in% 'meaningful']
  nonclust.tbl = pos.counts[cluster %in% 'nonMeaningful']

  if(nrow(clust.tbl) == 0){
    #message(paste('No meaningful positions found for', gene, sep=' '))
    return(NULL)
  }

  clust.tbl$distance = c(0,diff(clust.tbl$pos)) #calculate inter event distance.

  #If more than one meaningful positions are found within a 5 aa distance, join them to form a cluster.
  if(nrow(clust.tbl) > 1){

    #initialize variables.
    cstart = end = clust.tbl[1,pos]
    n = clust.tbl[1,N]
    cdf = c()
    cluster = 1

    #Go through entire table and update variables.
    for(i in 2:nrow(clust.tbl)){
      pos = clust.tbl[i,pos]

      d = clust.tbl[i,distance]

      if(d < mergeDist){
        end = pos
        n = n + clust.tbl[i,N]
      }else{
        tempdf = data.frame(cluster = paste('cluster', cluster, sep='_'), start = cstart, end = end ,N = n)
        cdf = rbind(cdf, tempdf)
        cstart = end = pos
        n = clust.tbl[i,N]
        cluster = cluster + 1
      }
    }
    cdf = rbind(cdf, data.frame(cluster = paste('cluster', cluster, sep='_'), start = cstart, end = end ,N = n))
  } else {
    cdf = data.frame(cluster = 'cluster_1', start = clust.tbl$pos, end = clust.tbl$pos ,N = clust.tbl$N)
  }

  #merge adjacent variants to clusters.
  for(i in 1:nrow(cdf)){
    tempcdf = cdf[i,]
    nonclust.tbl$startDist = nonclust.tbl$pos - tempcdf$start
    nonclust.tbl$endDist = nonclust.tbl$pos - tempcdf$end

    merge.adj.to.start = nonclust.tbl[startDist >= -5 & startDist <= 0]
    if(nrow(merge.adj.to.start) > 0){
      tempcdf$start = merge.adj.to.start[which(merge.adj.to.start$startDist == min(merge.adj.to.start$startDist)),pos]
      tempcdf$N = tempcdf$N + sum(merge.adj.to.start$N)
    }

    merge.adj.to.end = nonclust.tbl[endDist <= 5 & endDist >= 0]
    if(nrow(merge.adj.to.end) > 0){
      tempcdf$end = merge.adj.to.end[which(merge.adj.to.end$endDist == max(merge.adj.to.end$endDist)),pos]
      tempcdf$N = tempcdf$N + sum(merge.adj.to.end$N)
    }
    cdf[i,] = tempcdf
  }
  cdf$Hugo_Symbol = gene

  #Calcluate cluster score.

  total.muts = nrow(prot.dat) #total variants for this gene.
  clusterScores = c()

  for(i in 1:nrow(cdf)){
    temp.prot.dat = prot.dat[pos >= as.numeric(cdf$start[i]) & pos <= as.numeric(cdf$end[i])]
    temp.prot.dat.summary = temp.prot.dat[,.N, pos]
    temp.prot.dat.summary[,fraction:= N/total.muts]

    peak = temp.prot.dat.summary[N == max(N), pos]

    posVector = as.numeric(temp.prot.dat.summary[,pos])
    fractionMutVector = unlist(lapply(posVector, FUN = function(x) temp.prot.dat.summary[pos == x, fraction]))
    distanceVector = suppressWarnings(abs(posVector - peak))

    clusterScores = c(clusterScores,  sum( fractionMutVector / (sqrt(2)^ distanceVector)))

  }

  cdf$clusterScore = clusterScores

  gene.clust.res = data.frame(Hugo_Symbol = gene, clusters = nrow(cdf), muts_in_clusters = sum(cdf$N), clusterScores = sum(cdf$clusterScore), protLen = protLen)
  return(gene.clust.res)
}

# -------------------End of clustering function-------------------
