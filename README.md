## maftools - An R package to summarize, analyze and visualize MAF files.

###Introduction. 

With advances in Cancer Genomics, Mutation Annotation Format (MAF) is being widley accepted and used to store variants detected. 
[The Cancer Genome Atlas](http://cancergenome.nih.gov) Project has seqenced over 30 different cancers with sample size of each cancer type being over 200. The [resulting data](https://wiki.nci.nih.gov/display/TCGA/TCGA+MAF+Files) consisting of genetic variants is stored in the form of [Mutation Annotation Format](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification). 
This package attempts to summarize, analyze, annotate and visualize MAF files in an efficient manner either from TCGA sources or any in-house studies as long as the data is in MAF format. Maftools can also handle ICGC Simple Somatic Mutation format.

maftools is on :point_right: [bioRxiv](http://biorxiv.org/content/early/2016/05/11/052662) :bowtie:

Please cite the below if you find this tool useful for you.

Mayakonda, A. and H.P. Koeffler, Maftools: Efficient analysis, visualization and summarization of MAF files from large-scale cohort based cancer studies. bioRxiv, 2016. doi: http://dx.doi.org/10.1101/052662


### MAF field requirements.
MAF files contains many fields ranging from chromosome names to cosmic annotations. However, most of the analysis in maftools uses following fields. Please stick to MAF specifications for better results.
  * Mandatoty fields: __Hugo_Symbol, Chromosome, Start_Position, End_position, Variant_Classification, Variant_Type and Tumor_Sample_Barcode__. 
  * Recommended optional fields: non MAF specific fields containing vaf and amino acid change information. 
Complete specififcation of MAF files can be found on [NCI TCGA page](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification).

NOTE: If you have variants stored as VCFs or as an MAF like tab seperated format, convert them to MAF using [vcf2maf/maf2maf](https://github.com/mskcc/vcf2maf). Merge MAFs from all samples into a single MAF before processing with maftools.


## Vignette and a case study.
A complete documentation of maftools using TCGA LAML<sup>1</sup> as a case study can be found [here](http://www.bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html).

<p align="center">
<img src="https://github.com/PoisonAlien/PoisonAlien.github.io/blob/master/images/maftools.gif">
</p>



### Stuffs maftools can do.
1. Analysis
  * Detect cancer driver genes based on positional clustering of variants<sup>2</sup>.
  * Detect Mutually exclusive set of genes<sup>3</sup>.
  * Compare two MAF files (cohorts) to detect differentially mutated genes.
  * Add pfam domains and summarize.
  * Extract mutational signatures and compare them to validated signatures<sup>4</sup>.
  * APOBEC Enrichment score estimation<sup>5</sup>.
  * Tumor heterogenity and MATH (Mutant-Allele Tumor Heterogeneity) score estimation<sup>6</sup>.
  * Read and summarize GISTIC results.
  * Pan-cancer analysis/comparisison
  * Survival analysis
  * Compare mutation load against all 33 TCGA cohorts
2. Rich Visualizations
  * Make oncoplots<sup>7</sup>.
  * Make lollipop plots.
  * Map variants on copy number (CBS) segments 
  * Forest plots
  * Plot Transitions and Transversions. 
  * Plot maf summary.
  * CoOncoplots
  * Genecloud
  * Rainfall plots and change point detection
3. Annotation
  * Annotate variants locally using Oncotator API.
  * Convert Annovar annotations into MAF.
  * Convert ICGC simple somatic mutation format into MAF.


#### Installation:

Easy way: Install from [Bioconductor](http://bioconductor.org/packages/release/bioc/html/maftools.html).

```{r}
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("maftools")
```

Install from Github for updated features (some of functions from here may not be available on Bioconductor release branch).

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

For full documentation please refer to [vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html).


###References.
1.	Cancer Genome Atlas Research, N., Genomic and epigenomic landscapes of adult de novo acute myeloid leukemia. N Engl J Med, 2013. 368(22): p. 2059-74.
2.	Tamborero, D., A. Gonzalez-Perez, and N. Lopez-Bigas, OncodriveCLUST: exploiting the positional clustering of somatic mutations to identify cancer genes. Bioinformatics, 2013. 29(18): p. 2238-44.
3.	Leiserson, M.D., Wu, H.T., Vandin, F. & Raphael, B.J. CoMEt: a statistical approach to identify combinations of mutually exclusive alterations in cancer. Genome Biol 16, 160 (2015).
4.	Alexandrov, L.B., et al., Signatures of mutational processes in human cancer. Nature, 2013. 500(7463): p. 415-21.
5. Roberts SA, Lawrence MS, Klimczak LJ, et al. An APOBEC Cytidine Deaminase Mutagenesis Pattern is Widespread in Human Cancers. Nature genetics. 2013;45(9):970-976. doi:10.1038/ng.2702.
6.	Mroz, E.A. & Rocco, J.W. MATH, a novel measure of intratumor genetic heterogeneity, is high in poor-outcome classes of head and neck squamous cell carcinoma. Oral Oncol 49, 211-5 (2013).
7.	Gu, Z., Eils, R. & Schlesner, M. Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics (2016).

#### Powered By:
* [data.table](https://github.com/Rdatatable/data.table/wiki) at [warp speed](https://en.wikipedia.org/wiki/Warp_drive)
* [ggplot2](https://github.com/hadley/ggplot2)
