<img src="vignettes/maftools_hex.svg" align="left" height="140" /></a>

## maftools - An R package to summarize, analyze and visualize MAF files

[![bioc](http://www.bioconductor.org/shields/downloads/maftools.svg)](https://bioconductor.org/packages/stats/bioc/maftools/) 
[![bioc](http://www.bioconductor.org/shields/years-in-bioc/maftools.svg)](http://bioconductor.org/packages/devel/bioc/html/maftools.html)
[![bioc](http://www.bioconductor.org/shields/build/devel/bioc/maftools.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/maftools/)
[![GitHub closed issues](https://img.shields.io/github/issues-closed-raw/poisonalien/maftools.svg)](https://github.com/poisonalien/maftools/issues)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Lines Of Code](https://tokei.rs/b1/github/poisonalien/maftools?category=code)](https://github.com/poisonalien/maftools)

## Introduction

`maftools` provides a comprehensive set of functions for processing [MAF](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) files and to perform most commonly used analyses in cancer genomics. See [here](http://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html) for a detailed usage and a case study.

## Installation
```{r}
#Install from Bioconductor repository
BiocManager::install("maftools")

#Install from GitHub repository
BiocManager::install("PoisonAlien/maftools")
```

## Getting started: Vignette and a case study

A complete documentation of maftools using [TCGA LAML](https://www.nejm.org/doi/full/10.1056/nejmoa1301689) as a case study can be found [here](http://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html).

<p align="left">
<img src="https://user-images.githubusercontent.com/8164062/97981605-d8a59500-1dd2-11eb-9f5e-cc808f7b3f91.gif" height="320">
</p>

## Citation

**_Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP. 2018. Maftools: efficient and comprehensive analysis of somatic variants in cancer. [Genome Research](https://doi.org/10.1101/gr.239244.118). PMID: [30341162](https://www.ncbi.nlm.nih.gov/pubmed/?term=30341162)_**


## Useful links


|                                                    File formats                                                    |                                           Data Portals                                          |                                        Annotation tools                                       |
|:------------------------------------------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|:---------------------------------------------------------------------------------------------:|
|               [Mutation Annotation Format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)              |                               [TCGA](http://cancergenome.nih.gov)                               |       [vcf2maf](https://github.com/mskcc/vcf2maf) - for converting your VCF files to MAF      |
|                      [Variant Call Format](https://en.wikipedia.org/wiki/Variant_Call_Format)                      |                                  [ICGC](https://docs.icgc.org/)                                 | Ensembl Variant Effect Predictor [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html) |
| ICGC [Simple Somatic Mutation Format](https://docs.icgc.org/submission/guide/icgc-simple-somatic-mutation-format/) |                        [Broad Firehose](https://gdac.broadinstitute.org/)                       |           [Annovar](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/)           |
|                                                                                                                    |                            [cBioPortal](https://www.cbioportal.org/)                            |    [Funcotator](https://gatk.broadinstitute.org/hc/en-us/articles/360037224432-Funcotator)    |
|                                                                                                                    |        [CIViC](https://civicdb.org/home) - Clinical interpretation of variants in cancer        |                                                                                               |
|                                                                                                                    | [DGIdb](http://www.dgidb.org/) - Information on drug-gene interactions and the druggable genome |                                                                                               |


## Similar packages/tools

Below are some more useful software packages for somatic variant analysis (not necessarily similar to maftools)

* [TRONCO](https://github.com/BIMIB-DISCo/TRONCO) - Repository of the TRanslational ONCOlogy library (R)
* [dndscv](https://github.com/im3sanger/dndscv) - dN/dS methods to quantify selection in cancer and somatic evolution (R)
* [cloneevol](https://github.com/hdng/clonevol) - Inferring and visualizing clonal evolution in multi-sample cancer sequencing (R)
* [sigminer](https://github.com/ShixiangWang/sigminer) - Primarily for signature analysis and visualization in R. Supports `maftools` output (R)
* [GenVisR](https://github.com/griffithlab/GenVisR) - Primarily for visualization (R)
* [comut](https://github.com/vanallenlab/comut) - Primarily for visualization (Python)
* [TCGAmutations](https://github.com/PoisonAlien/TCGAmutations) - pre-compiled curated somatic mutations from TCGA cohorts (from Broad Firehose and TCGA MC3 Project) that can be loaded into `maftools` (R)

***

#### Powered By

* [data.table](https://github.com/Rdatatable/data.table/wiki) at [warp speed](https://en.wikipedia.org/wiki/Warp_drive)
