getCounts=function(var,bam, MapQuality, BaseQuality, ref_genome){


  tot.df=read.delim('var_temp.txt', header=F, stringsAsFactors = F, sep='\t', as.is = T, colClasses = c(rep('character',5)))
  rownames(tot.df)=paste(tot.df$V1,tot.df$V2,sep=":")
  #tot.df = tot.df[,1:5]

  system(command = paste("bam-readcount -q",MapQuality," -i -w1 -b",BaseQuality," -l var_temp.txt -f ",ref_genome,bam,"> dx.counts"), ignore.stdout = F, ignore.stderr = F)
  dx.counts=read.delim(pipe("cut -f 1-10 dx.counts"),header=F)
  rownames(dx.counts)=paste(dx.counts$V1, dx.counts$V2,sep = ":")

  tot.df = tot.df[rownames(dx.counts),]

  tot.df$depth=dx.counts$V4
  tot.df$A=sapply(strsplit(as.character(dx.counts$V6),split = ":",fixed = T),"[",2)
  tot.df$C=sapply(strsplit(as.character(dx.counts$V7),split = ":",fixed = T),"[",2)
  tot.df$G=sapply(strsplit(as.character(dx.counts$V8),split = ":",fixed = T),"[",2)
  tot.df$T=sapply(strsplit(as.character(dx.counts$V9),split = ":",fixed = T),"[",2)
  tot.df$N=sapply(strsplit(as.character(dx.counts$V10),split = ":",fixed = T),"[",2)

  dx.vaf=c()
  dx.count=c()

  for(i in 1:nrow(tot.df)) {
    alt=as.character(tot.df[i,5])
    if(alt == "A"){
      dx.vaf=c(dx.vaf,as.numeric(as.character(tot.df[i,7]))/as.numeric(as.character(tot.df[i,6])))
      dx.count=c(dx.count,as.numeric(as.character(tot.df[i,7])))
    }

    if(alt == "C"){
      dx.vaf=c(dx.vaf,as.numeric(as.character(tot.df[i,8]))/as.numeric(as.character(tot.df[i,6])))
      dx.count=c(dx.count,as.numeric(as.character(tot.df[i,8])))
    }

    if(alt == "G"){
      dx.vaf=c(dx.vaf,as.numeric(as.character(tot.df[i,9]))/as.numeric(as.character(tot.df[i,6])))
      dx.count=c(dx.count,as.numeric(as.character(tot.df[i,9])))
    }

    if(alt == "T"){
      dx.vaf=c(dx.vaf,as.numeric(as.character(tot.df[i,10]))/as.numeric(as.character(tot.df[i,6])))
      dx.count=c(dx.count,as.numeric(as.character(tot.df[i,10])))
    }
  }

  tot.df$vaf=dx.vaf*100
  colnames(tot.df)[1:5]=c("chr","start","end","ref","alt")

  dx.ref.count=c()

  for(i in 1:nrow(tot.df)) {
    ref=as.character(tot.df[i,4])
    if(ref == "A"){
      dx.ref.count=c(dx.ref.count,as.numeric(as.character(tot.df[i,7])))

    } else if(ref == "C"){
      dx.ref.count=c(dx.ref.count,as.numeric(as.character(tot.df[i,8])))

    } else if(ref == "G"){
      dx.ref.count=c(dx.ref.count,as.numeric(as.character(tot.df[i,9])))

    } else if(ref == "T"){
      dx.ref.count=c(dx.ref.count,as.numeric(as.character(tot.df[i,10])))
    }
  }

  prim.tum=cbind(tot.df[,c(1:3,6)],dx.ref.count,dx.count,tot.df[,12])
  colnames(prim.tum)[4:7]=c("t_depth","t_ref_count","t_alt_count","t_vaf")
  #prim.tum$chr = gsub(pattern = "chr",replacement = "",x = prim.tum$chr)
  system(command = 'rm dx.counts')
  system(command = 'rm var_temp.txt')
  return(prim.tum[,c(1,2,4:7)])
}
