sampleSwaps = function(bams = NULL, build = "hg19", prefix = NULL, add = TRUE, min_depth = 30){

  build = match.arg(arg = build, choices = c("hg19", "hg38"))

  if(build == "hg19"){
    snps = system.file("extdata", "hg19_smash_snps.tsv.gz", package = "maftools")
    #snps = "inst/extdata/hg19_smash_snps.tsv.gz"
  }else{
    #snps = "inst/extdata/hg38_smash_snps.tsv.gz"
    snps = system.file("extdata", "hg38_smash_snps.tsv.gz", package = "maftools")
  }
  snps = data.table::fread(input = snps, sep = "\t")

  if(!is.null(prefix)){
    if(add){
      snps$chr = paste(prefix, snps$chr, sep = '')
    }else{
      snps$chr = gsub(pattern = prefix, replacement = '', x = snps$chr, fixed = TRUE)
    }
  }

  rc = bamreadcounts(bam = bams, loci = snps)

  snps[, id := paste0(chr, ":", start)]

  rc_af = lapply(rc, function(x){
    xrc = merge(x, snps[,.(id, ref, alt)], by.x = 'loci', by.y = 'id')

    xrc_acounts = apply(xrc, 1, function(x){

      ref_allele = x['ref']
      ref_rc = 0

      ref_rc = switch (ref_allele,
                       A = x['A'],
                       T = x['T'],
                       G = x['G'],
                       C = x['C'])

      alt_allele = x['alt']
      alt_rc = 0

      alt_rc = switch (alt_allele,
                       A = x['A'],
                       T = x['T'],
                       G = x['G'],
                       C = x['C'])

      vaf_tbl = data.table::data.table(ref_rc = as.numeric(ref_rc), alt_rc = as.numeric(alt_rc), loci = x["loci"])
      vaf_tbl[, vaf := alt_rc/(ref_rc+alt_rc)]
      vaf_tbl
    })
    xrc = merge(xrc, data.table::rbindlist(l = xrc_acounts), by = "loci")
    xrc[,.(loci, ref_rc, alt_rc, vaf)]
  })

  rc_bind = data.table::rbindlist(l = rc_af, idcol = "sample")
  rc_bind = rc_bind[!is.nan(vaf)]
  rc_bind = rc_bind[, total := ref_rc + alt_rc][total > min_depth]

  rc_df = data.table::dcast(data = rc_bind, loci ~ sample, value.var = 'vaf', fill = NA)
  data.table::setDF(x = rc_df, rownames = rc_df$loci)
  rc_df$loci = NULL
  rc_df = rc_df[complete.cases(rc_df),]

  list(rc_df, rc_af)
}
