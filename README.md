### maftools - An R package to summarize, analyze and visualize MAF files.

####Introduction.
With advances in Cancer Genomics, maf format is being widley accepted and used to store variants detected. 
[The Cancer Genome Atlas](http://cancergenome.nih.gov) Project has seqenced over 30 different cancers with sample size of each cancer type being over 200. The [resulting data](https://wiki.nci.nih.gov/display/TCGA/TCGA+MAF+Files) consisting of genetic variants is stored in the form of [Mutation Annotation Format](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification). 
This package attempts to summarize, analyze, annotate and visualize MAF files in an efficient manner either from TCGA sources or any in-house studies as long as the data is in MAF format.

#### MAF field requirements.
MAF files contains many fields ranging from chromosome names to cosmic annotations. However, most of the analysis in maftools uses following fields.
  * Mandatoty fields: __Hugo_Symbol, Chromosome, Start_Position, End_position, Variant_Classification, Variant_Type and Tumor_Sample_Barcode__. 
  * Recommended optional fields: non MAF specific fields containing vaf and amino acid change information. 
Complete specififcation of MAF files can be found on [NCI TCGA page](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification).

#### Vignette and a case study.
A complete documentation of maftools using TCGA LAML<sup>1</sup> as a case study can be found [here](http://poisonalien.github.io).

#### Stuffs maftools can do.
1. Analysis
  * Detect cancer driver genes based on positional clustering of variants<sup>2</sup>.
  * Detect Mutually exclusive set of genes<sup>5</sup>.
  * Extract mutational signatures and compare them to validated signatures.
  * Tumor heterogenity and MATH (Mutant-Allele Tumor Heterogeneity) score estimation.
  * Add pfam domains and summarize.
2. Visualization
  * Make oncoplots.
  * Make lollipop plots.
  * Plot Transitions and Transversions. 
  * Plot maf summary.
3. Annotation
  * Annotate variants locally using Oncotator API.
  * Convert Annovar annotations into MAF.

#### Installation:

```{r results='hide'}
#Install Bioconductor dependencies.
source("http://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
biocLite("VariantAnnotation")
biocLite("Biostrings")

#Install maftools from github repository.
library("devtools")
install_github(repo = "PoisonAlien/maftools")
```

#### Reading MAF files
First we read maf file using function `read.maf` which reads maf, summarises variants in various ways, creates oncomatrix and performs sorting.

```{r results='hide'}
require(maftools)
#read TCGA maf file for LAML
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
laml = read.maf(maf = laml.maf, removeSilent = T, useAll = F)
```
MAF file is read and stored as an MAF object of S4 class. Each slot can be accessed using directly `@`. However there are accessor methods which can be used also.

```{r, echo=TRUE, size=1}
#MAF file summarized by Genes 
getGeneSummary(laml)

#MFA file summarized by samples (Tumor Sample Barcode)
getSampleSummary(laml)

#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = laml, basename = 'laml')

```
#### Quicky plot maf stats

```{r, echo=TRUE}
plotmafSummary(maf = laml, rmOutlier = T, addStat = 'median')
```
![image1](https://github.com/PoisonAlien/maftoolsDump/tree/master/images/image1.png)

#### oncoplot to summarize maf file
`oncoplot` This function uses modified oncoprint script from [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) package by [Zuguang Gu](https://github.com/jokergoo) with added functionalities, while taking care of format conversions.

```{r, echo=TRUE, fig.height=7,fig.width=12}
#We will plot top ten mutated genes
oncoplot(maf = laml, top = 10)
```
![image2](https://github.com/PoisonAlien/maftoolsDump/tree/master/images/image2)

Use arguments drawRowBar, drawColBar to control side and upper barplots.

#### Adding annotations to oncoplot
We can add annotations to the bottom of the plot.
```{r, echo=TRUE, fig.height=9,fig.width=12}
#Read FAB classification of TCGA LAML barodes.
laml.fab.anno = system.file('extdata', 'tcga_laml_fab_annotation.txt', package = 'maftools')
laml.fab.anno = read.delim(laml.fab.anno, sep = '\t')
head(laml.fab.anno)
#We will plot same top ten mutated genes with FAB classification as annotation.
oncoplot(maf = laml, top = 10, annotation = laml.fab.anno)
```
![image3](https://github.com/PoisonAlien/maftoolsDump/tree/master/images/image3)

#### Classify SNVs into Trasitions and Transversions
Each Single Nucleotide Variant can be classified into [Trasition or Transversion]((http://www.mun.ca/biology/scarr/Transitions_vs_Transversions.html)). Variants can also be divided into six different classes, which helps to know us which kind of conversions are more frequent in a given type of cancer.  

```{r, echo=TRUE,fig.height=4,fig.width=6, warning=FALSE,fig.align='center'}
laml.titv = titv(maf = laml, useSyn = T)
```
![image5](https://github.com/PoisonAlien/maftoolsDump/tree/master/images/image5.png)

It also returns a list of dataframes with raw counts for each conversion, fraction of each conversion and Ti to Tv ratios.

##Lollipop plots for amino acid changes.
We can map protein changes on to the Protein structure similar to those draw by [ProteinPaint](https://pecan.stjude.org/proteinpaint/TP53/) or [MutationMapper](http://www.cbioportal.org/mutation_mapper.jsp). Labelling requires [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html) package.
```{r, echo = TRUE, fig.height=4,fig.width=7,fig.align='center'}
lollipopPlot(maf = laml, gene = 'KIT', AACol = 'Protein_Change', labelPos = 'all')
```
![image8](https://github.com/PoisonAlien/maftoolsDump/tree/master/images/image8)

####Detcting cancer causing genes.
maftools comes with the function `oncodrive` which identifies cancer genes (driver) from a given MAF. `oncodrive` is a based on algorithm [oncodriveCLUST](http://bg.upf.edu/group/projects/oncodrive-clust.php) which was originally implemented in Python. Concept is based on the fact that most of the variants in cancer causing genes are enriched at few specific loci (aka hotspots). This method takes advantage of such positions to identify cancer genes. If you use this function, please cite [OncodriveCLUST article](http://bioinformatics.oxfordjournals.org/content/early/2013/07/31/bioinformatics.btt395.full)<sup>2</sup>.

```{r}
laml.sig = oncodrive(maf = aml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = T)
```
![image10](https://github.com/PoisonAlien/maftoolsDump/tree/master/images/image10)

####Summarizing by pfam domains.
`pfamDomain` adds domain information and summarizes amino acid changes accoriding to the domains that are affected. This serves the puposes of knowing what domain in given cancer cohort, is most frequently affected. This function is inspired from Pfam annotation modulce of MuSic tool<sup>3</sup>.

```{r}
laml.pfam = pfamDomains(maf = laml, AACol = 'Protein_Change', top = 10)
#Protein summary (Printing first 7 columns for display convenience)
laml.pfam$proteinSummary[,1:7, with = F]
#Domain summary (Printing first 3 columns for display convenience)
laml.pfam$domainSummary[,1:3, with = F]
```
![image11](https://github.com/PoisonAlien/maftoolsDump/tree/master/images/image11.png)

#### Annotating variants with Oncotator
We can also annotate variants using [oncotator](http://www.broadinstitute.org/oncotator/) API<sup>4</sup>.

```{r}
var.file = system.file('extdata', 'variants.tsv', package = 'maftools')
#This is what input looks like
var = read.delim(var.file, sep = '\t')
head(var)
```

```{r, results='hide'}
#Annotate 
var.maf = oncotate(maflite = var.file, header = T)
```

This is quite time consuming if input is big.

#### Mutual Exclusivity and Oncoprint.
Many genes in cancer show strong exclusiveness in mutation pattern. We can detect such pair of genes using `mutExclusive` which runs `comet_exact_test` from `cometExactTest` package for significance<sup>5</sup>. 

```{r, echo = TRUE, fig.height=1.5,fig.width=7,fig.align='center'}
mutExclusive(maf = laml, genes = c('NPM1', 'RUNX1'))
#Output
##   n.00 n.01 n.10 n.11 gene1 gene2                pval
## 1  124   18   54    0  NPM1 RUNX1 0.00109823009964897
#We can visualize this pair using oncostrip. `oncostrip` draws a matrix similar to [oncoprint](http://www.cbioportal.org/faq.jsp#what-are-oncoprints) on [cBioPortal](http://www.cbioportal.org/index.do).

oncostrip(maf = laml, genes = c('NPM1', 'RUNX1'), sort = T, legend = T, removeNonMutated = T)
```
![image6](https://github.com/PoisonAlien/maftoolsDump/tree/master/images/image6)

####Tumor Heterogenity
Tumors are generally heterogenous i.e, consist of multiple clones. This heterogenity can be inferred by clustering variant allele frequencies. We will manually mention vaf column. Requires [mclust](https://cran.r-project.org/web/packages/mclust/index.html) package. Although mlcust performs fairly well, there are other tools like [SciClone](https://github.com/genome/sciclone) which does better job at clustering and density estimation<sup>5</sup>.

```{r, echo = TRUE, fig.align='center', fig.height=5, fig.width=7}
#We will run this for sample TCGA.AB.2972
inferHetrogentiy(maf = laml, tsb = 'TCGA.AB.2972', vafCol = 'TumorVAF_WU')
```
![image7](https://github.com/PoisonAlien/maftoolsDump/tree/master/images/image7)

#### Extract Mutation Signatures
Every cancer, as it progresses leaves a signature characterised by specific pattern of nucleotide substitutions. [Alexandrov et.al](http://www.nature.com/nature/journal/v500/n7463/full/nature12477.html) have shown such signatures, derived from over 7000 cancer samples. Such signatures can be extracted by decomposiong matrix of nucleotide substitutions, classified into 96 substitution classes based on immediate bases sorrouding the mutated base. Extracted signatures can also be compared to those [21 validated signatures](http://cancer.sanger.ac.uk/cosmic/signatures)<sup>6</sup>. 

ExtractSignatures uses `nmf` from [NMF](https://cran.r-project.org/web/packages/NMF/index.html) package to decompose the matrix and extract signatures<sup>7</sup>. 
NOTE: Reading fasta file is memory consuming and it occupies ~3gb of memory while extracting adjacent bases from human genome.

Trying this on TCGA Liver MAF file.
```{r}
#Read MAF file
maf = "hgsc.bcm.edu_LIHC.IlluminaGA_DNASeq.1.somatic.maf"
lihc = read.maf(maf = maf, removeSilent = T, useAll = F)
#Extract adjacent bases from immediate 3' and 5', classify them into 96 conversion classes.
lihc.tnm = trinucleotideMatrix(maf = lihc, ref_genome = 'hg19.fa', prefix = 'chr', add = T, useSyn = T)
## reading fasta (this might take a while)..
## Extracting adjacent bases..
## matrix of dimension 198x96

# Decompose extracted matrix using Non-negative Matrix Factorization.
lihc.signatures = extractSignatures(mat = lihc.tnm)
## Estimating best rank..
## Using 2 as a best-fit rank based on maximum cophenetic correlation coefficient.
## Comparing against experimentally validated 21 signatures.. (See Alexandrov et.al Nature 2013 for details.)
## Found Signature_1 most similar to validated Signature_4. Correlation coeff: 0.547099730968907 
## Found Signature_2 most similar to validated Signature_12. Correlation coeff: 0.647694192108207 
```
Signature_2 which corelates will validated Signature_12 was observed in Liver samples characterised by T>C mutations showing transcriptional strand bias<sup>8</sup>.

![image9](https://github.com/PoisonAlien/maftoolsDump/tree/master/images/image9)

#### Add read count and allele frequencies to maf.
`addReadCounts()` adds read depths for reference and alternate allele from corresponding bam file. This internally runs [bam-readcount](https://github.com/genome/bam-readcount) to get the counts and adds them to maf file. 

#### Other functions
For full documentation please refer to [vignette](http://poisonalien.github.io).

####References.
1.	Cancer Genome Atlas Research, N., Genomic and epigenomic landscapes of adult de novo acute myeloid leukemia. N Engl J Med, 2013. 368(22): p. 2059-74.
2.	Tamborero, D., A. Gonzalez-Perez, and N. Lopez-Bigas, OncodriveCLUST: exploiting the positional clustering of somatic mutations to identify cancer genes. Bioinformatics, 2013. 29(18): p. 2238-44.
3.	Dees, N.D., et al., MuSiC: identifying mutational significance in cancer genomes. Genome Res, 2012. 22(8): p. 1589-98.
4.	Ramos, A.H., et al., Oncotator: cancer variant annotation tool. Hum Mutat, 2015. 36(4): p. E2423-9.
5.	Leiserson, M.D., et al., CoMEt: a statistical approach to identify combinations of mutually exclusive alterations in cancer. Genome Biol, 2015. 16: p. 160.
6.	Miller, C.A., et al., SciClone: inferring clonal architecture and tracking the spatial and temporal patterns of tumor evolution. PLoS Comput Biol, 2014. 10(8): p. e1003665.
7.	Alexandrov, L.B., et al., Signatures of mutational processes in human cancer. Nature, 2013. 500(7463): p. 415-21.
8.	Gaujoux, R. and C. Seoighe, A flexible R package for nonnegative matrix factorization. BMC Bioinformatics, 2010. 11: p. 367.
