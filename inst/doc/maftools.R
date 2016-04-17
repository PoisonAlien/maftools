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
laml = read.maf(maf = laml.maf, removeSilent = TRUE, useAll = FALSE)

## ------------------------------------------------------------------------
#Typing laml shows basic summary of MAF file.
laml
#Shows sample summry.
getSampleSummary(laml)
#Shows frequently mutated genes.
getGeneSummary(laml)
#Shows all fields in MAF
getFields(laml)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = laml, basename = 'laml')

## ----fig.height=5, fig.width=7-------------------------------------------
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median')

## ---- fig.align='left',fig.height=6,fig.width=10, eval=T, fig.align='left'----
#We will draw oncoplots for top ten mutated genes. (Removing non-mutated samples from the plot for better visualization)
oncoplot(maf = laml, top = 10, removeNonMutated = TRUE)

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
oncoplot(maf = laml, top = 10, annotation = laml.fab.anno, removeNonMutated = TRUE, colors = col)

## ---- fig.height=2,fig.width=8,fig.align='center'------------------------
oncostrip(maf = laml, genes = c('DNMT3A','NPM1', 'RUNX1'), removeNonMutated = TRUE, showTumorSampleBarcodes = FALSE)

## ---- fig.align='default', fig.height=6, fig.width=8, eval = T-----------
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)

## ----fig.height=4,fig.width=7,fig.align='center'-------------------------
#Lets plot lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
dnmt3a.lpop = lollipopPlot(maf = laml, gene = 'DNMT3A', AACol = 'Protein_Change')

## ----fig.height=3,fig.width=8,fig.align='center'-------------------------
#Lets mutations on KIT gene, without repel option.
kit.lpop = lollipopPlot(maf = laml, gene = 'KIT', AACol = 'Protein_Change', labelPos = c(416, 418), refSeqID = 'NM_000222')
#Same plot with repel=T
kit.lpop = lollipopPlot(maf = laml, gene = 'KIT', AACol = 'Protein_Change', labelPos = c(416, 418), refSeqID = 'NM_000222', repel = TRUE)

## ---- fig.height=4,fig.width=8,fig.align='center'------------------------
tcga.ab.009.seg <- system.file("extdata", "TCGA.AB.3009.hg19.seg.txt", package = "maftools")
plotCBSsegments(cbsFile = tcga.ab.009.seg, maf = laml, labelAll = TRUE)

## ------------------------------------------------------------------------
#We will run mutExclusive on top 10 mutated genes. 
laml.mut.excl = mutExclusive(maf = laml, top = 10)
head(laml.mut.excl)

## ---- fig.height=1.5,fig.width=8,fig.align='center'----------------------
oncostrip(maf = laml, genes = c('NPM1', 'RUNX1'), sort = TRUE, removeNonMutated = TRUE)

## ---- fig.align='default', fig.width=7,fig.height=5, message=F,results='hide', eval=T----
laml.sig = oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
head(laml.sig)

## ---- fig.align='default', fig.width=7,fig.height=5, eval= T-------------
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE)

## ---- fig.align='left',fig.width=7, fig.height=5, eval=T-----------------
laml.pfam = pfamDomains(maf = laml, AACol = 'Protein_Change', top = 10)
#Protein summary (Printing first 7 columns for display convenience)
laml.pfam$proteinSummary[,1:7, with = FALSE]
#Domain summary (Printing first 3 columns for display convenience)
laml.pfam$domainSummary[,1:3, with = FALSE]

## ----results='hide', message=FALSE---------------------------------------
#Primary APL MAF
primary.apl = system.file("extdata", "APL_primary.maf.gz", package = "maftools")
primary.apl = read.maf(maf = primary.apl)
#Relapse APL MAF
relapse.apl = system.file("extdata", "APL_relapse.maf.gz", package = "maftools")
relapse.apl = read.maf(maf = relapse.apl)

## ------------------------------------------------------------------------
#We will consider only genes which are mutated in at-least in 5 samples in one of the cohort, to avoid single mutated genes.
pt.vs.rt <- mafCompare(m1 = primary.apl, m2 = relapse.apl, m1Name = 'PrimaryAPL', m2Name = 'RelapseAPL', minMut = 5)
print(pt.vs.rt)

## ---- fig.width=6, fig.height=6, fig.align='center'----------------------
apl.pt.vs.rt.fp = forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, show = 'stat', color = c('royalblue', 'maroon'))

## ---- echo = TRUE, fig.align='center', fig.height=5, fig.width=7, eval=T----
#We will run this for sample TCGA.AB.2972
tcga.ab.2972.het = inferHeterogeneity(maf = laml, tsb = 'TCGA.AB.2972', vafCol = 'i_TumorVAF_WU')
print(tcga.ab.2972.het$clusterMeans)
#Visualizing results
plotClusters(clusters = tcga.ab.2972.het)

## ---- fig.align='center', fig.height=5, fig.width=7, eval=T--------------
seg = system.file('extdata', 'TCGA.AB.3009.hg19.seg.txt', package = 'maftools')
tcga.ab.3009.het = inferHeterogeneity(maf = laml, tsb = 'TCGA.AB.3009', segFile = seg)
#Visualizing results. Highlighting those variants on copynumber altered variants.
plotClusters(clusters = tcga.ab.3009.het, genes = 'CN_altered', showCNvars = TRUE)

## ---- results='hide', message=F, fig.align='center',fig.width=7, fig.height=6, eval = T----
#we will specify for random 4 patients.
laml.math = math.score(maf = laml, vafCol = 'i_TumorVAF_WU', 
                       sampleName = c('TCGA.AB.3009', 'TCGA.AB.2849', 'TCGA.AB.3002', 'TCGA.AB.2972'))

## ---- eval=T-------------------------------------------------------------
print(laml.math)

## ---- eval=F-------------------------------------------------------------
#  #First we extract adjacent bases to the mutated locus and clssify them into 96 substitution classes.
#  laml.tnm = trinucleotideMatrix(maf = laml, ref_genome = '/path/to/hg19.fa',
#                                 prefix = 'chr', add = TRUE, ignoreChr = 'chr23', useSyn = TRUE)

## ---- fig.height=5, fig.width=5, eval=F----------------------------------
#  #Run main function with maximum 6 signatures.
#  laml.sign = extractSignatures(mat = laml.tnm, nTry = 6, plotBestFitRes = FALSE)
#  ## Warning : Found zero mutations for conversions A[T>G]C
#  ## Estimating best rank..
#  ## Using 2 as a best-fit rank based on maximum cophenetic correlation coefficient.
#  ## Comparing against experimentally validated 21 signatures.. (See Alexandrov et.al Nature 2013 for details.)
#  ## Found Signature_1 most similar to validated Signature_1B. Correlation coeff: 0.788340817468816
#  ## Found Signature_2 most similar to validated Signature_1A. Correlation coeff: 0.724579532825333

## ---- echo=F-------------------------------------------------------------
laml.sign = structure(list(signatures = structure(c(0.0188253312203968, 0.0108560245997858, 
0.00570899206006542, 9.02584597832905e-19, 9.02584597832905e-19, 
9.02584597832905e-19, 0.00798908345795425, 9.02584597832905e-19, 
9.02584597832905e-19, 9.02584597832905e-19, 0.00554013015871254, 
0.00906586191997992, 0.0129424152140228, 0.00915581596661287, 
9.02584597832905e-19, 0.0294145800318701, 0.0152955816165724, 
9.02584597832905e-19, 0.00705949920764882, 9.02584597832905e-19, 
9.02584597832905e-19, 0.00250097818128397, 0.017648748019122, 
9.02584597832905e-19, 0.0117085058837841, 0.00705949920764882, 
0.00145423657396816, 0.00705949920764882, 0.00351261465262355, 
0.00823608240892362, 0.00655428772997421, 3.16584451150982e-10, 
0.0447101616484425, 0.0168776517741237, 0.0843006539138647, 9.02584597832905e-19, 
0.0507204981749981, 0.037243056940506, 0.096122157993921, 0.0405193900343099, 
0.0447101616484425, 0.0200288705914744, 0.027962792487326, 0.0333515555926063, 
0.022203427900867, 9.02584597832905e-19, 0.0495594817896151, 
0.0233849028874635, 0.00214023037309487, 0.00210330402362136, 
9.02584597832905e-19, 9.02584597832905e-19, 0.00235316640254961, 
0.00823608240892362, 0.0141189984152976, 0.00651077698107633, 
0.00206331531764248, 0.0105892488114732, 0.00588291600637401, 
9.02584597832905e-19, 0.0011765832012748, 9.02584597832905e-19, 
0.000647524758589879, 0.00470633280509921, 0.00800218554138202, 
9.02584597832905e-19, 0.0118522019873281, 0.0329443296356945, 
0.00331719478032513, 0.00780753137981535, 9.02584597832905e-19, 
0.0211784976229465, 9.02584597832905e-19, 0.00468206699897619, 
0.00707759003249691, 0.00474683447685859, 0.00235316640254961, 
0.00470633280509921, 0.00636590445820346, 9.02584597832905e-19, 
0.00470633280509921, 0.0100038703594569, 9.02584597832905e-19, 
9.02584597832905e-19, 9.02584597832905e-19, 9.02584597832905e-19, 
0.00514936857903451, 0.00588291600637401, 0.00352974960382441, 
9.02584597832905e-19, 9.02584597832905e-19, 9.02584597832905e-19, 
0.00352974960382441, 9.02584597832905e-19, 0.00235316640254961, 
9.02584597832905e-19, 8.08587587136962e-19, 0.00913575679424068, 
0.00387726069897548, 0.00747816360942334, 0.0121520158653129, 
0.0130867863164908, 0.00580485814591106, 0.0121520158653129, 
0.0149563272188467, 0.00747816360942334, 0.00120710655212577, 
0.00962323270360145, 8.08587587136962e-19, 0.0123560774709514, 
0.011217245414135, 8.08587587136962e-19, 8.08587587136962e-19, 
0.00560862270706751, 8.08587587136962e-19, 0.00747816360942334, 
0.0130867863164909, 0.00268687802374693, 8.08587587136962e-19, 
0.00841293406060126, 0.0019150852985174, 8.08587587136962e-19, 
0.00445326252778329, 8.08587587136962e-19, 0.00188315425842584, 
8.08587587136962e-19, 0.004140461630987, 0.00934770426025953, 
8.08587587136962e-19, 0.039872978595969, 0.0723057119036414, 
0.0233692612794479, 0.0242027988241808, 0.0685620768790476, 0.0965655133550832, 
0.0883935904074667, 8.08587587136962e-19, 0.0289564678801311, 
0.0778045914738293, 0.0108937130102545, 0.0188158940614483, 0.0289778839865154, 
0.0297990540254919, 0.0160076962880105, 0.00203871407099345, 
0.00861144443235405, 0.00280431135353375, 0.00280431135353375, 
8.08587587136962e-19, 8.08587587136962e-19, 8.08587587136962e-19, 
0.00137071864655311, 0.00209982145465293, 8.08587587136962e-19, 
8.08587587136962e-19, 0.00467385225588959, 8.08587587136962e-19, 
0.00186954090235584, 0.00322463708879178, 8.08587587136962e-19, 
0.0188812351344225, 0.0130867863164908, 0.00927908539498494, 
8.08587587136962e-19, 0.00577749301564168, 0.00781863833182143, 
0.0177606385723804, 8.08587587136962e-19, 0.00841293406060126, 
0.00375836047409369, 0.00652902038079363, 0.00370690408371161, 
8.08587587136962e-19, 8.08587587136962e-19, 0.00242058725094486, 
0.0186954090235583, 8.08587587136962e-19, 0.00953234416709177, 
0.00186954090235583, 0.00280431135353375, 0.000934770451177918, 
0.00841293406060126, 0.00525664016147213, 8.08587587136962e-19, 
8.08587587136962e-19, 0.00280431135353375, 0.00467385225588959, 
0.00373908180471167, 8.08587587136962e-19, 0.00467385225588959, 
8.08587587136962e-19, 0.00560862270706751), .Dim = c(96L, 2L), .Dimnames = list(
    c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", 
    "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", 
    "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "A[C>G]A", 
    "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", 
    "C[C>G]T", "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", 
    "T[C>G]C", "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", 
    "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", 
    "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", 
    "T[C>T]T", "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", 
    "C[T>A]C", "C[T>A]G", "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", 
    "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "A[T>C]A", 
    "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", 
    "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", 
    "T[T>C]C", "T[T>C]G", "T[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", 
    "A[T>G]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", 
    "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", 
    "T[T>G]T"), c("Signature_1", "Signature_2"))), corTable = structure(c(0.778439049788447, 
0.724579532825333, 0.788340817468816, 0.718314020391453, 0.0774615504710396, 
0.061187479472836, -0.0405758863668983, -0.113253571934054, 0.0215437598720942, 
0.0368280494167392, 0.397886860567507, 0.321437635997587, 0.628887363450396, 
0.649246431686657, 0.231410563208561, 0.375463335348244, 0.198815299931041, 
0.133061340212547, -0.00839716809342338, -0.0688476751301835, 
0.253697010273925, 0.0654371338622379, 0.281411060546338, 0.511188585647839, 
-0.0126454799211963, 0.0273671159443648, 0.00350075761425314, 
0.00611361803498541, 0.369554253203226, 0.430341741741454, 0.281135513302693, 
0.423803590022172, 0.0865853537467677, 0.0184651013556118, 0.0234382417598754, 
-0.0288268638451162, 0.0524982989107897, 0.0340916078412234, 
0.551305287258253, 0.653122978837916, 0.350607320032417, 0.306339259742027, 
0.0812573289808097, 0.0855544290480513), .Dim = c(2L, 22L), .Dimnames = list(
    c("Signature_1", "Signature_2"), c("Signature_1A", "Signature_1B", 
    "Signature_2", "Signature_3", "Signature_4", "Signature_5", 
    "Signature_6", "Signature_7", "Signature_8", "Signature_9", 
    "Signature_10", "Signature_11", "Signature_12", "Signature_13", 
    "Signature_14", "Signature_15", "Signature_16", "Signature_17", 
    "Signature_18", "Signature_19", "Signature_20", "Signature_21"
    )))), .Names = c("signatures", "corTable"))

## ---- fig.width=7, fig.height=5, fig.align='center', eval = T------------
plotSignatures(laml.sign)

## ------------------------------------------------------------------------
var.file = system.file('extdata', 'variants.tsv', package = 'maftools')
#This is what input looks like
var = read.delim(var.file, sep = '\t')
head(var)

## ---- results='hide', eval=F, message=F----------------------------------
#  #Annotate
#  var.maf = oncotate(maflite = var.file, header = TRUE)

## ---- eval = F-----------------------------------------------------------
#  #Results from oncotate. First 20 columns.
#  var.maf[1:10, 1:20, with = FALSE]

## ---- eval=T-------------------------------------------------------------
var.annovar = system.file("extdata", "variants.hg19_multianno.txt", package = "maftools")
var.annovar.maf = annovarToMaf(annovar = var.annovar, Center = 'CSI-NUS', refBuild = 'hg19', 
                               tsbCol = 'Tumor_Sample_Barcode', table = 'ensGene', header = TRUE)

print(var.annovar.maf)

## ------------------------------------------------------------------------
##Extract data for samples 'TCGA.AB.3009' and 'TCGA.AB.2933'  (Printing just 5 rows for display convenience)
subsetMaf(maf = laml, tsb = c('TCGA.AB.3009', 'TCGA.AB.2933'))[1:5]
##Same as above but return output as an MAF object
subsetMaf(maf = laml, tsb = c('TCGA.AB.3009', 'TCGA.AB.2933'), mafObj = TRUE)

## ------------------------------------------------------------------------
##Select all Splice_Site mutations from DNMT3A and NPM1
subsetMaf(maf = laml, genes = c('DNMT3A', 'NPM1'), query = "Variant_Classification == 'Splice_Site'")
##Select all variants with VAF above 30%
subsetMaf(maf = laml, query = "t_vaf > 30")
##Same as above but include only 'Protein_Change' column in the output
subsetMaf(maf = laml, query = "t_vaf > 30", fields = 'AAChange')

