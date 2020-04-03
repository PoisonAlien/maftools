<img src="https://github.com/PoisonAlien/maftools/blob/master/vignettes/maftools.png" height="90" width="240" />

## maftools - An R package to summarize, analyze and visualize MAF files

[![bioc](http://www.bioconductor.org/shields/downloads/maftools.svg)](https://bioconductor.org/packages/stats/bioc/maftools/) 
[![bioc](http://www.bioconductor.org/shields/years-in-bioc/maftools.svg)](http://bioconductor.org/packages/devel/bioc/html/maftools.html)
[![bioc](http://www.bioconductor.org/shields/build/devel/bioc/maftools.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/maftools/)
[![GitHub issues](https://img.shields.io/github/issues-raw/poisonalien/maftools.svg)](https://github.com/poisonalien/maftools/issues)
[![GitHub closed issues](https://img.shields.io/github/issues-closed-raw/poisonalien/maftools.svg)](https://github.com/poisonalien/maftools/issues)
[![Github Stars](https://img.shields.io/github/stars/poisonalien/maftools.svg?style=social&label=Github)](https://github.com/poisonalien/maftools)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)


## Introduction
With advances in Cancer Genomics, Mutation Annotation Format (MAF) is being widley accepted and used to store variants detected. 
[The Cancer Genome Atlas](http://cancergenome.nih.gov) Project has seqenced over 30 different cancers with sample size of each cancer type being over 200. The [resulting data](https://wiki.nci.nih.gov/display/TCGA/TCGA+MAF+Files) consisting of genetic variants is stored in the form of [Mutation Annotation Format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/). 
This package attempts to summarize, analyze, annotate and visualize MAF files in an efficient manner either from TCGA sources or any in-house studies as long as the data is in MAF format. Maftools can also handle ICGC Simple Somatic Mutation format.

## Citation

***
**_Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP. 2018. Maftools: efficient and comprehensive analysis of somatic variants in cancer. [Genome Research](https://doi.org/10.1101/gr.239244.118). PMID: [30341162](https://www.ncbi.nlm.nih.gov/pubmed/?term=30341162)_**

***

### Installation

```{r}
#Install from Bioconductor repository
BiocManager::install("maftools")

#Install from GitHub repository
BiocManager::install("PoisonAlien/maftools")
```

## Vignette and a case study
A complete documentation of maftools using TCGA LAML<sup>1</sup> as a case study can be found [here](http://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html).

<p align="center">
<img src="https://github.com/PoisonAlien/PoisonAlien.github.io/blob/master/images/maftools.gif">
</p>


## Usage
Simple usage: Just read maf file using `read.maf` and pass the resulting maf object to any one of the function for plotting or analysis or set operations.

<p align="center">
<img src="https://github.com/PoisonAlien/maftools/blob/master/vignettes/overview.png">
</p>

### Stuffs maftools can do
1. Analysis
  * Detect cancer driver genes based on positional clustering of variants<sup>2</sup>
  * Detect Mutually exclusive set of genes<sup>3</sup>
  * Compare two MAF files (cohorts) to detect differentially mutated genes
  * Add pfam domains and summarize
  * Extract mutational signatures and compare them to validated signatures<sup>4</sup>
  * APOBEC Enrichment score estimation<sup>5</sup>
  * Tumor heterogenity and MATH (Mutant-Allele Tumor Heterogeneity) score estimation<sup>6</sup>
  * Read and summarize GISTIC results
  * Pan-cancer analysis/comparisison
  * Survival analysis
  * Compare mutation load against all 33 TCGA cohorts
2. Rich Visualizations
  * Make oncoplots
  * Make lollipop plots
  * Map variants on copy number (CBS) segments 
  * Forest plots
  * Plot Transitions and Transversions
  * Plot maf summary
  * CoOncoplots
  * Genecloud
  * Rainfall plots and change point detection
3. Annotation
  * Annotate variants locally using Oncotator API
  * Convert Annovar annotations into MAF
  * Convert ICGC simple somatic mutation format into MAF

### References

1.	Cancer Genome Atlas Research, N., Genomic and epigenomic landscapes of adult de novo acute myeloid leukemia. N Engl J Med, 2013. 368(22): p. 2059-74.
2.	Tamborero, D., A. Gonzalez-Perez, and N. Lopez-Bigas, OncodriveCLUST: exploiting the positional clustering of somatic mutations to identify cancer genes. Bioinformatics, 2013. 29(18): p. 2238-44.
3.	Leiserson, M.D., Wu, H.T., Vandin, F. & Raphael, B.J. CoMEt: a statistical approach to identify combinations of mutually exclusive alterations in cancer. Genome Biol 16, 160 (2015).
4.	Alexandrov, L.B., et al., Signatures of mutational processes in human cancer. Nature, 2013. 500(7463): p. 415-21.
5. Roberts SA, Lawrence MS, Klimczak LJ, et al. An APOBEC Cytidine Deaminase Mutagenesis Pattern is Widespread in Human Cancers. Nature genetics. 2013;45(9):970-976. doi:10.1038/ng.2702.
6.	Mroz, E.A. & Rocco, J.W. MATH, a novel measure of intratumor genetic heterogeneity, is high in poor-outcome classes of head and neck squamous cell carcinoma. Oral Oncol 49, 211-5 (2013).

#### Powered By
* [data.table](https://github.com/Rdatatable/data.table/wiki) at [warp speed](https://en.wikipedia.org/wiki/Warp_drive)
