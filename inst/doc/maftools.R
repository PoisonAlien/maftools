## ----results='hide', eval=F----------------------------------------------
#  #Install Bioconductor dependencies.
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("ComplexHeatmap")
#  biocLite("VariantAnnotation")
#  biocLite("Biostrings")
#  
#  #Install maftools from github repository.
#  library("devtools")
#  install_github(repo = "PoisonAlien/maftools")

## ----results='hide', message=FALSE---------------------------------------
suppressWarnings(require(maftools))
#read TCGA maf file for LAML
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
laml = read.maf(maf = laml.maf, removeSilent = T, useAll = F)

## ------------------------------------------------------------------------
#Typing laml shows basic summary of MAF file.
laml
#Shows sample summry.
getSampleSummary(laml)
#Shows frequently mutated genes.
getGeneSummary(laml)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = laml, basename = 'laml')

## ----fig.height=5, fig.width=7-------------------------------------------
plotmafSummary(maf = laml, rmOutlier = T, addStat = 'median')

## ---- fig.align='left',fig.height=6,fig.width=10, eval=T, fig.align='left'----
#We will draw oncoplots for top ten mutated genes. (Removing non-mutated samples from the plot for better visualization)
oncoplot(maf = laml, top = 10, removeNonMutated = T)

## ---- fig.height=7,fig.width=10, eval=T, fig.align='left'----------------
#Read FAB classification of TCGA LAML barcodes.
laml.fab.anno = system.file('extdata', 'tcga_laml_fab_annotation.txt', package = 'maftools')
laml.fab.anno = read.delim(laml.fab.anno, sep = '\t')
head(laml.fab.anno)
#Changing colors (You can use any colors, here in this example we will use a color palette from RColorBrewer)
col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
#We will plot same top ten mutated genes with FAB classification as annotation and using above defined colors.
oncoplot(maf = laml, top = 10, annotation = laml.fab.anno, removeNonMutated = T, colors = col)

## ---- fig.height=2,fig.width=8,fig.align='center'------------------------
oncostrip(maf = laml, genes = c('DNMT3A','NPM1', 'RUNX1'), removeNonMutated = T, showTumorSampleBarcodes = F)

## ---- fig.align='default', fig.height=6, fig.width=8, eval = T-----------
laml.titv = titv(maf = laml, plot = F, useSyn = T)
#plot titv summary
plotTiTv(res = laml.titv)

## ----fig.height=4,fig.width=7,fig.align='center'-------------------------
#Lets plot lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
dnmt3a.lpop = lollipopPlot(maf = laml, gene = 'DNMT3A', AACol = 'Protein_Change')

## ----fig.height=4,fig.width=7,fig.align='center'-------------------------
#Lets mutations on KIT gene, without repel option.
kit.lpop = lollipopPlot(maf = laml, gene = 'KIT', AACol = 'Protein_Change', labelPos = c(416, 418), refSeqID = 'NM_000222')
#Same plot with repel=T
kit.lpop = lollipopPlot(maf = laml, gene = 'KIT', AACol = 'Protein_Change', labelPos = c(416, 418), refSeqID = 'NM_000222', repel = T)

## ------------------------------------------------------------------------
#We will run mutExclusive on top 10 mutated genes. 
laml.mut.excl = mutExclusive(maf = laml, top = 10)
head(laml.mut.excl)

## ---- fig.height=1.5,fig.width=8,fig.align='center'----------------------
oncostrip(maf = laml, genes = c('NPM1', 'RUNX1'), sort = T, removeNonMutated = T)

## ---- fig.align='default', fig.width=7,fig.height=5, message=F,results='hide', eval=T----
laml.sig = oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')

## ---- fig.align='default', fig.width=7,fig.height=5, eval= T-------------
head(laml.sig)
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = T)

## ---- fig.align='left',fig.width=7, fig.height=5, eval=T-----------------
laml.pfam = pfamDomains(maf = laml, AACol = 'Protein_Change', top = 10)
#Protein summary (Printing first 7 columns for display convenience)
laml.pfam$proteinSummary[,1:7, with = F]
#Domain summary (Printing first 3 columns for display convenience)
laml.pfam$domainSummary[,1:3, with = F]

## ---- echo = TRUE, fig.align='center', fig.height=5, fig.width=7, eval=T----
#We will run this for sample TCGA.AB.2972
inferHeterogeneity(maf = laml, tsb = 'TCGA.AB.2972', vafCol = 'i_TumorVAF_WU')

## ---- results='hide', message=F, fig.align='center',fig.width=7, fig.height=6, eval = T----
#we will specify for random 4 patients.
laml.math = math.score(maf = laml, vafCol = 'i_TumorVAF_WU', 
                       sampleName = c('TCGA.AB.3009', 'TCGA.AB.2849', 'TCGA.AB.3002', 'TCGA.AB.2972'))

## ---- eval=T-------------------------------------------------------------
print(laml.math)

## ---- eval=T-------------------------------------------------------------
#First we extract adjacent bases to the mutated locus and clssify them into 96 substitution classes.
laml.tnm = trinucleotideMatrix(maf = laml, ref_genome = '~/NGS/gatk_ref/hg19.fa', 
                               prefix = 'chr', add = T, ignoreChr = 'chr23', useSyn = T)

## ---- fig.height=5, fig.width=5, eval=T----------------------------------
#Run main function with maximum 6 signatures. 
laml.sign = extractSignatures(mat = laml.tnm, nTry = 6)

## ---- fig.width=7, fig.height=5, fig.align='center', eval = T------------
plotSignatures(laml.sign)

## ------------------------------------------------------------------------
var.file = system.file('extdata', 'variants.tsv', package = 'maftools')
#This is what input looks like
var = read.delim(var.file, sep = '\t')
head(var)

## ---- results='hide', eval=F, message=F----------------------------------
#  #Annotate
#  var.maf = oncotate(maflite = var.file, header = T)

## ---- eval = F-----------------------------------------------------------
#  #Results from oncotate. First 20 columns.
#  var.maf[1:10, 1:20, with =F]

## ---- eval=T-------------------------------------------------------------
var.annovar = system.file("extdata", "variants.hg19_multianno.txt", package = "maftools")
var.annovar.maf = annovarToMaf(annovar = var.annovar, Center = 'CSI-NUS', refBuild = 'hg19', 
                               tsbCol = 'Tumor_Sample_Barcode', table = 'ensGene', header = T)

print(var.annovar.maf)

