## The following script creates the ../inst/extdata/tcga_cohort.txt.gz dataset
## Requires TCGAmutations R package (https://github.com/PoisonAlien/TCGAmutations) and Agilent sureselect target_regions BED file (https://kb.10xgenomics.com/hc/en-us/articles/115004150923-Where-can-I-find-the-Agilent-Target-BED-files-)

#Get TCGA 33 cohort names
tcga_cohorts = TCGAmutations::tcga_available()[,Study_Abbreviation][1:33]
cohort_tissue_sites = c(
  ACC = "Adrenal Gland Carcinoma",
  BLCA = "Bladder",
  BRCA = "Breast",
  CESC = "Cervix",
  CHOL = "Bile Duct",
  COAD = "Colorectal Colon",
  DLBC = "Lymph Nodes",
  ESCA = "Esophagus",
  GBM = "Brain Glioblastoma",
  HNSC = "Head and Neck",
  KICH = "Kidney Chromophobe",
  KIRC = "Kidney Clear Cell",
  KIRP = "Kidney Papilary",
  LAML = "Bone Marrow",
  LGG = "Brain Glioma",
  LIHC = "Liver",
  LUAD = "Lung Adeno",
  LUSC = "Lung Squamous",
  MESO = "Pleura",
  OV = "Ovary",
  PAAD = "Pancreas",
  PCPG = "Adrenal Gland",
  PRAD = "Prostate",
  READ = "Colorectal Rectum",
  SARC = "Soft Tissue",
  SKCM = "Skin",
  STAD = "Stomach",
  TGCT = "Testis",
  THCA = "Thyroid",
  THYM = "Thymus",
  UCEC = "Uterus Corpus",
  UCS = "Uterus",
  UVM = "Eye"
)

#Get Agilent target regions (Has to be downloaded from Agilent site, requires registration)
#I downloaded it on Oct 6, 2020, for hg19.
ss50 = data.table::fread(input = "~/Dropbox/Agilent_SS_Human_All_Exon_V7/S31285117_Regions.bed", skip = 2)
ss50$V1 = gsub(pattern = "chr", replacement = "", x = ss50$V1)

#target region size (should be 35804808 i.e, OR 35.8 MB)
sum(ss50[,V3 - V2])/1e6

tcga_tmb = lapply(tcga_cohorts, function(cohort) {
  cohort_maf = TCGAmutations::tcga_load(study = cohort, source = "MC3")
  cohort_maf_ss50 = subsetMaf(maf = cohort_maf, ranges = ss50)
  cohort_tmb = getSampleSummary(cohort_maf_ss50)[, .(Tumor_Sample_Barcode, total)]
  cohort_tmb$site = cohort_tissue_sites[cohort]
  cohort_tmb
})

names(tcga_tmb) = tcga_cohorts
tcga_tmb = data.table::rbindlist(l = tcga_tmb, use.names = TRUE, fill = TRUE, idcol = "cohort")
data.table::fwrite(x = tcga_tmb, file = "inst/extdata/tcga_cohort.txt.gz", sep = "\t")
