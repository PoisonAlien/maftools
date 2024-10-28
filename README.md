<img src="vignettes/maftools_hex.svg" align="left" height="140" /></a>

## maftools - An R package to summarize, analyze and visualize MAF files

[![GitHub closed issues](https://img.shields.io/github/issues-closed-raw/poisonalien/maftools.svg)](https://github.com/poisonalien/maftools/issues)
[![R-CMD-check](https://github.com/PoisonAlien/maftools/workflows/R-CMD-check/badge.svg)](https://github.com/PoisonAlien/maftools/actions)

## Introduction

maftools is a comprehensive toolkit for processing somatic variants from cohort-based cancer genomic studies. maftools offers over 80 functions to perform the most commonly required tasks in cancer genomics, using [MAF](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) as the only input file type.

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
<img src="https://user-images.githubusercontent.com/8164062/97981605-d8a59500-1dd2-11eb-9f5e-cc808f7b3f91.gif" height="320" height="400">
</p>

## Primary applications 

maftools is extremely easy to use, starting with importing an [MAF](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) file along with the associated clinical data. Once the data is successfully imported, the resulting MAF object can be passed to various functions. Key applications include:

- [Cohort summarization using oncoplots](https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/oncoplots.html#08_Combining_everything)
- [Identify co-occurring and mutually exclusive events](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#91_Somatic_Interactions)
- [Clinical enrichment analysis](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#96_Clinical_enrichment_analysis)
- [Detect cancer driver genes](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#92_Detecting_cancer_driver_genes_based_on_positional_clustering)
- [Infer tumor heterogeneity](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#99_Tumor_heterogeneity_and_MATH_scores)
- [Analyze known cancer signaling pathways](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#98_Oncogenic_Signaling_Pathways)
- [De-novo somatic signature analysis with NMF](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#9103_Signature_analysis)
- [Compare two cohorts to identify differentially mutated genes](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#95_Comparing_two_cohorts_(MAFs))
- [Perform survival analysis and predict genesets associated with survival](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#942_Predict_genesets_associated_with_survival)
- [Drug-gene interactions](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#97_Drug-Gene_Interactions)

Besides the MAF files, maftools can handle sequencing alignment BAM files, copy number output from GISTIC and mosdepth. Please refer to the package documentation sections below to learn more.

- [Generate personalized cancer report](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/cancer_hotspots.html) for known somatic [hotspots](https://www.cancerhotspots.org/)
- [Sample mismatch and relatedness analysis](https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html#12_Sample_swap_identification)
- [Copy number analysis](https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/cnv_analysis.html) with [ASCAT](https://github.com/VanLoo-lab/ascat) and [mosdepth](https://github.com/brentp/mosdepth)

Moreover, analyzing all 33 TCGA cohorts along with the harmonized clinical data is a breeze. 

- A single command [tcgaLoad](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#13_TCGA_cohorts) will import the desired TCGA cohort thereby avoiding costly time spent on data mining from public databases. 
- Please refer to an associated software package [TCGAmutations](https://github.com/PoisonAlien/TCGAmutations) that provides ready to use `MAF` objects for 33 TCGA cohorts and 2427 cell line profiles from CCLE - along with relevant clinical information for all sequenced samples.

## Citation

**_Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP. 2018. Maftools: efficient and comprehensive analysis of somatic variants in cancer. [Genome Research](https://doi.org/10.1101/gr.239244.118). PMID: [30341162](https://www.ncbi.nlm.nih.gov/pubmed/?term=30341162)_**


## Useful links

| File Fomats                                                                                                        | Data portals                                                                                    | Annotation tools                                                                                                                       |
|--------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------|
| [Mutation Annotation Format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)                            | [TCGA](http://cancergenome.nih.gov)                                                             | [vcf2maf](https://github.com/mskcc/vcf2maf) - for converting your VCF files to MAF                                                     |
| [Variant Call Format](https://en.wikipedia.org/wiki/Variant_Call_Format)                                           | [ICGC](https://docs.icgc.org/)                                                                  | [annovar2maf](https://github.com/PoisonAlien/annovar2maf) - for converting annovar output files to MAF                                 |
| ICGC [Simple Somatic Mutation Format](https://docs.icgc.org/submission/guide/icgc-simple-somatic-mutation-format/) | [Broad Firehose](https://gdac.broadinstitute.org/)                                              | [bcftools csq](https://samtools.github.io/bcftools/howtos/csq-calling.html) - Rapid annotations of VCF files with variant consequences |
|                                                                                                                    | [cBioPortal](https://www.cbioportal.org/)                                                       | [Annovar](https://annovar.openbioinformatics.org/en/latest/)                                                              |
|                                                                                                                    | [PeCan](https://pecan.stjude.cloud/)                                                            | [Funcotator](https://gatk.broadinstitute.org/hc/en-us/articles/360037224432-Funcotator)                                                |
|                                                                                                                    | [CIViC](https://civicdb.org/home) - Clinical interpretation of variants in cancer               |                                                                                                                                        |
|                                                                                                                    | [DGIdb](http://www.dgidb.org/) - Information on drug-gene interactions and the druggable genome |                                                                                                                                        |


## Useful packages/tools

Below are some more useful software packages for somatic variant analysis

* [TRONCO](https://github.com/BIMIB-DISCo/TRONCO) - Repository of the TRanslational ONCOlogy library (R)
* [dndscv](https://github.com/im3sanger/dndscv) - dN/dS methods to quantify selection in cancer and somatic evolution (R)
* [cloneevol](https://github.com/hdng/clonevol) - Inferring and visualizing clonal evolution in multi-sample cancer sequencing (R)
* [sigminer](https://github.com/ShixiangWang/sigminer) - Primarily for signature analysis and visualization in R. Supports `maftools` output (R)
* [GenVisR](https://github.com/griffithlab/GenVisR) - Primarily for visualization (R)
* [comut](https://github.com/vanallenlab/comut) - Primarily for visualization (Python)
* [TCGAmutations](https://github.com/PoisonAlien/TCGAmutations) - pre-compiled curated somatic mutations from TCGA cohorts (from Broad Firehose and TCGA MC3 Project) that can be loaded into `maftools` (R)
* [somaticfreq](<https://github.com/PoisonAlien/somaticfreq>) - rapid genotyping of known somatic hotspot variants from the tumor BAM files. Generates a browsable/sharable HTML report. (C)

***

#### Powered By

* [data.table](https://github.com/Rdatatable/data.table/wiki) at [warp speed](https://en.wikipedia.org/wiki/Warp_drive)
