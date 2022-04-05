#Code to create SNP6 loci used by `gtSNPs()`
#Files are located at: https://github.com/CompEpigen/ezASCAT/blob/main/inst/extdata/GRCh*_SNP6.tsv.gz

library(data.table)

#hg38
download.file(url = "https://api.gdc.cancer.gov/data/77fbfff6-2acc-47ca-a5f6-c488beb46879", destfile = "snp6.na35.liftoverhg38.txt.zip")
unzip(zipfile = "snp6.na35.liftoverhg38.txt.zip")
snp6_hg38 = data.table::fread(input = "snp6.na35.liftoverhg38.txt")
snp6_hg38 = snp6_hg38[!type %in% "CN"][order(as.character(chr), as.numeric(pos))][,.(chr, pos, probeid)]
nrow(snp6_hg38) #932,148 loci

#hg19: Get `GPL6801-4019.txt` from GEO annotations https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6801
## On build 32 (ucsc=hg19; ncbi=GRCh37)
snp6_hg19 = data.table::fread(input = "GPL6801-4019.txt", skip = 10)
snp6_hg19 = snp6_hg19[grepl(pattern = "^SNP", x = snp6_hg19$ID)][!Chromosome %in% "---"][order(Chromosome, as.numeric(`Physical Position`))]
nrow(snp6_hg19) #930,104 loci

length(intersect(x = snp6_hg19$ID, snp6_hg38$probeid)) #be 929,132 probes common

data.table::fwrite(x = snp6_hg19[,.(Chromosome, `Physical Position`)], file = "inst/extdata/GRCh37_SNP6.tsv.gz", sep = "\t", col.names = FALSE)
data.table::fwrite(x = snp6_hg38[,.(chr, pos)], file = "inst/extdata/GRCh38_SNP6.tsv.gz", sep = "\t", col.names = FALSE)

system(command = "rm snp6.na35.liftoverhg38.txt")
system(command = "rm snp6.na35.liftoverhg38.txt.zip")
system(command = "rm GPL6801-4019.txt")
