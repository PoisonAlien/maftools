# SNPs for hg19 and hg38 are compiled from the publication:
#
# Westphal, M., Frankhouser, D., Sonzone, C. et al.
# SMaSH: Sample matching using SNPs in humans. BMC Genomics 20, 1001 (2019).
# https://doi.org/10.1186/s12864-019-6332-7

hg19_snps = data.table::fread(input = "https://raw.githubusercontent.com/rbundschuh/SMaSH/master/snps_hg19.vcf", sep = "\t")[,.(V1, V2, V4, V5, V3)]
colnames(hg19_snps) = c('chr', 'start', 'ref', 'alt', 'rsid')

hg38_snps = data.table::fread(input = "https://raw.githubusercontent.com/rbundschuh/SMaSH/master/snps_GRCh38.vcf", sep = "\t")[,.(V1, V2, V4, V5, V3)]
colnames(hg38_snps) = c('chr', 'start', 'ref', 'alt', 'rsid')

data.table::fwrite(x = hg19_snps, file = "inst/extdata/hg19_smash_snps.tsv.gz", sep = "\t", compress = "gzip")
data.table::fwrite(x = hg38_snps, file = "inst/extdata/hg38_smash_snps.tsv.gz", sep = "\t", compress = "gzip")
