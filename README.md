### maftools - A package to process and visualize maf files. 

With advances in Cancer Genomics, maf format is being widley accepted and used to store variants detected. 
[The Cancer Genome Atlas](http://cancergenome.nih.gov) Project has seqenced over 30 different cancers with sample size of each cancer type being over 200. The [resulting data](https://wiki.nci.nih.gov/display/TCGA/TCGA+MAF+Files) consisting of genetic variants is stored in the form of [Mutation Annotation Format](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification). This package attempts to summerize such files either from TCGA or in house studies by providing various functions for plotting and parsing.

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

#### Reading MAF file
First we read maf file using fuction `read.maf` which also summarises variants by various ways and sorts them. 

```{r results='hide'}
require(maftools)
#read TCGA maf file for LAML
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
laml = read.maf(maf = laml.maf, removeSilent = T, useAll = F)
```

MAF file is read and stored as an MAF object of S4 class. Each slot can be accessed using directly @. However there are accessor methods which can be used also.


```{r, echo=TRUE, size=1}
#MAF file summarized by Genes 
getGeneSummary(laml)

#MFA file summarized by samples (Tumor Sample Barcode)
getSampleSummary(laml)

```
#### Quicky plot maf stats

```{r, echo=TRUE}
plotmafSummary(laml)
```
![image1](https://github.com/PoisonAlien/maftools/blob/master/images/image1)

#### oncoplot to summarize maf file
`oncoplot` This function uses slightly modified [oncoprint](https://github.com/jokergoo/ComplexHeatmap/blob/908b32ee4c495c74adfa077c967024a77c56b375/vignettes/oncoprint.R) script from [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) package by [Zuguang Gu](https://github.com/jokergoo), while taking care of format conversions and with some added fucntionalities.

```{r, echo=TRUE, fig.height=7,fig.width=14}
#We will plot top ten mutated genes
oncoplot(maf = laml, top = 10)
```
![image2](https://github.com/PoisonAlien/maftools/blob/master/images/image2)

Use arguments drawRowBar, drawColBar to control side and upper barplots.

#### Adding annotations to oncoplot
We can add annotations to the bottom of the plot.
```{r, echo=TRUE, fig.height=9,fig.width=16}
#Read FAB classification of TCGA LAML barodes.
laml.fab.anno = system.file('extdata', 'tcga_laml_fab_annotation.txt', package = 'maftools')
laml.fab.anno = read.delim(laml.fab.anno, sep = '\t')
head(laml.fab.anno)
#We will plot same top ten mutated genes with FAB classification as annotation.
oncoplot(maf = laml, top = 10, annotation = laml.fab.anno)
```
![image3](https://github.com/PoisonAlien/maftools/blob/master/images/image3)

#### Classify SNVs into Trasitions and Transversions
Each Single Nucleotide Variant can be classified into [Trasition or Transversion]((http://www.mun.ca/biology/scarr/Transitions_vs_Transversions.html)). Variants can also be divided into six different classes, which helps to know us which kind of conversions are more frequent in a given type of cancer.  

```{r, echo=TRUE,fig.height=4,fig.width=6, warning=FALSE,fig.align='center'}
laml.titv.summary = titv(maf = laml, plot = T)
```
![image5](https://github.com/PoisonAlien/maftools/blob/master/images/image5)

It also returns a list of dataframes with raw counts for each conversion, fraction of each conversion and Ti to Tv ratios.

##Lollipop plots for amino acid changes.
We can map protein changes on to the Protein structure similar to those draw by [ProteinPaint](https://pecan.stjude.org/proteinpaint/TP53/) or [MutationMapper](http://www.cbioportal.org/mutation_mapper.jsp). Labelling requires [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html) package.
```{r, echo = TRUE, fig.height=4,fig.width=7,fig.align='center'}
lollipopPlot(maf = laml, gene = 'KIT', label = T)
```
![image8](https://github.com/PoisonAlien/maftools/blob/master/images/image8)

####Detcting cancer causing genes.
maftools comes with the function `oncodrive` which identifies cancer genes (driver) from a given MAF. `oncodrive` is a based on algorithm [oncodriveCLUST](http://bg.upf.edu/group/projects/oncodrive-clust.php) which was originally implemented in Python. Concept is based on the fact that most of the variants in cancer causing genes are enriched at few specific loci (aka hotspots). This method takes advantage of such positions to identify cancer genes. If you use this function, please cite [OncodriveCLUST article](http://bioinformatics.oxfordjournals.org/content/early/2013/07/31/bioinformatics.btt395.full).

```{r}
laml.sig = oncodrive(maf = aml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = T)
```
![image10](https://github.com/PoisonAlien/maftools/blob/master/images/image10.png)

#### Annotating variants with Oncotator
We can also annotate variants using [oncotator](http://www.broadinstitute.org/oncotator/) API.

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
Many genes in cancer show strong exclusiveness in mutation pattern. We can detect such pair of genes using `mutExclusive` which runs `comet_exact_test` from `cometExactTest` package for significance. 

```{r, echo = TRUE, fig.height=1.5,fig.width=7,fig.align='center'}
mutExclusive(maf = laml, genes = c('NPM1', 'RUNX1'))

#Output
##   n.00 n.01 n.10 n.11 gene1 gene2                pval
## 1  124   18   54    0  NPM1 RUNX1 0.00109823009964897

#We can visualize this pair using oncoprint. `oncoprint` draws a matrix similar to [oncoprint](http://www.cbioportal.org/faq.jsp#what-are-oncoprints) on [cBioPortal](http://www.cbioportal.org/index.do).

oncoprint(maf = laml, genes = c('NPM1', 'RUNX1'), sort = T, legend = T, removeNonMutated = T)
```
![image6](https://github.com/PoisonAlien/maftools/blob/master/images/image6)

####Tumor Heterogenity
Tumors are generally heterogenous i.e, consist of multiple clones. This heterogenity can be inferred by clustering variant allele frequencies. We will manually mention vaf column. Requires [mclust](https://cran.r-project.org/web/packages/mclust/index.html) package. Although mlcust performs fairly well, it is recommended to try [SciClone](https://github.com/genome/sciclone) which is far superior for clustering and density estimation.

```{r, echo = TRUE, fig.align='center', fig.height=5, fig.width=7}
#We will run this for sample TCGA.AB.2972
inferHetrogentiy(maf = laml, tsb = 'TCGA.AB.2972', vafCol = 'TumorVAF_WU')
```
![image7](https://github.com/PoisonAlien/maftools/blob/master/images/image7)

#### Extract Mutation Signatures
Every cancer, as it progresses leaves a signature characterised by specific pattern of nucleotide substitutions. [Alexandrov et.al](http://www.nature.com/nature/journal/v500/n7463/full/nature12477.html) have shown such signatures, derived from over 7000 cancer samples. Such signatures can be extracted by decomposiong matrix of nucleotide substitutions, classified into 96 substitution classes based on immediate bases sorrouding the mutated base. Extracted signatures can also be compared to those [21 validated signatures](http://cancer.sanger.ac.uk/cosmic/signatures). 

ExtractSignatures uses nmf from [NMF](https://cran.r-project.org/web/packages/NMF/index.html) package to extract signatures. 
NOTE: Reading fasta file is memory consuming and it occupies ~3gb of memory while extracting adjacent bases from human genome.

Trying this on TCGA Liver MAF file.
```{r}
#Read MAF file
maf = "hgsc.bcm.edu_LIHC.IlluminaGA_DNASeq.1.somatic.maf"
lihc = read.maf(maf = maf, removeSilent = T, useAll = F)
#Extract adjacent bases from immediate 3' and 5', classify them into 96 conversion classes.
lihc.tnm = trinucleotideMatrix(maf = lihc, ref_genome = 'hg19.fa', prefix = 'chr', add = T)
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
Signature_2 which corelates will validated Signature_12 was observed in Liver samples characterised by T>C mutations showing transcriptional strand bias.

![image9](https://github.com/PoisonAlien/maftools/blob/master/images/image9)

#### Add read count and allele frequencies to maf.
`addReadCounts()` adds read depths for reference and alternate allele from corresponding bam file. This internally runs [bam-readcount](https://github.com/genome/bam-readcount) to get the counts and adds them to maf file. 
