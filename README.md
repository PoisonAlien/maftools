### maftools - A package to process and visualize maf files. 

With advances in Cancer Genomics, maf format is being widley accepted and used to store variants detected. 
[The Cancer Genome Atlas](http://cancergenome.nih.gov) Project has seqenced over 30 different cancers with sample size of each cancer type being over 200. The [resulting data](https://wiki.nci.nih.gov/display/TCGA/TCGA+MAF+Files) consisting of genetic variants is stored in the form of [Mutation Annotation Format](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification). This package attempts to summerize such files either from TCGA or in house studies by providing various functions for plotting and parsing.

#### Installation:

`library("devtools")`

`install_github(repo = "PoisonAlien/maftools")`

#### Dependencies: 
data.table, ggplot2, plyr, reshape, [cometExactTest](https://cran.r-project.org/web/packages/cometExactTest/)

Bioconductor packages:  [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html), [VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html) and [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html).

#### Reading MAF file
First we read maf file using fuction `read.maf` which also summarises variants by various ways and sorts them. It returns a list of `data.frame`s which can be accessed easily.

```{r results='hide'}
require(maftools)
#read TCGA maf file for LAML
laml.maf = system.file('extdata', 'tcga_laml.maf', package = 'maftools')
laml = read.maf(maf = laml.maf, removeSilent = T, useAll = F)
```

MAF file is summarized in various ways which can easily accessed. For example..

```{r, echo=TRUE, size=1}
#Based on Variant_Classification
laml$variant.classification.summary

#Based on frequcy of mutated geens
laml$gene.summary
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

#### oncoprint 
`oncoprint` draws a matrix similar to [oncoprint](http://www.cbioportal.org/faq.jsp#what-are-oncoprints) on [cBioPortal](http://www.cbioportal.org/index.do).

```{r, echo=TRUE,fig.height=1.5,fig.width=7,fig.align='center'}
oncoprint(maf = laml, genes = c('DNMT3A', 'NPM1'), sort = T, legend = T, removeNonMutated = T)
```
![image4](https://github.com/PoisonAlien/maftools/blob/master/images/image4)

#### Classify SNVs into Trasitions and Transversions
Each Single Nucleotide Variant can be classified into [Trasition or Transversion]((http://www.mun.ca/biology/scarr/Transitions_vs_Transversions.html)). Variants can also be divided into six different classes, which helps to know us which kind of conversions are more frequent in a given type of cancer.  

```{r, echo=TRUE,fig.height=4,fig.width=6, warning=FALSE,fig.align='center'}
laml.titv.summary = titv(maf = laml, plot = T)
```
![image5](https://github.com/PoisonAlien/maftools/blob/master/images/image5)

It also returns a list of dataframes with raw counts for each conversion, fraction of each conversion and Ti to Tv ratios.

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

#### Mutual Exclusivity
Many genes in cancer show strong exclusiveness in mutation pattern. We can detect such pair of genes using `mutExclusive` which runs `comet_exact_test` from `cometExactTest` package for significance. 

```{r, echo = TRUE, fig.height=1.5,fig.width=7,fig.align='center'}
mutExclusive(maf = laml, genes = c('NPM1', 'RUNX1'))
#We can visualize this pair using oncoprint
oncoprint(maf = laml, genes = c('NPM1', 'RUNX1'), sort = T, legend = T, removeNonMutated = T)
```
![image6](https://github.com/PoisonAlien/maftools/blob/master/images/image6)

#### Extract adjacent bases
One can also extract n number of adjacent (3' and 5') bases to the mutated locus using `addBases`. This is helpful in looking for [somatic-signatures](http://cancer.sanger.ac.uk/cosmic/signatures). This requires faidx indexed reference genome (fasta file).

#### Add read count and allele frequencies to maf.
`addReadCounts()` adds read depths for reference and alternate allele from corresponding bam file. This internally runs [bam-readcount](https://github.com/genome/bam-readcount) to get the counts and adds them to maf file. 
