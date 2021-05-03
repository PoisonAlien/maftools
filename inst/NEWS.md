# CHANGES IN VERSION 2.7.30 (GitHub/Bioc devel version)
## BUG FIXES
- `coOncoplot` not allowing more than one additional feature. Issue: [675](https://github.com/PoisonAlien/maftools/issues/675)

## ENHANCEMENTS
- Added protein domains for the gene `ALMS1`. Issue: [705](https://github.com/PoisonAlien/maftools/issues/705)
- Added `titv_col` argumtn to oncoplot. Issue: [702](https://github.com/PoisonAlien/maftools/issues/702)
- Added protein domains for the gene `FAM205A`. Issue: [701](https://github.com/PoisonAlien/maftools/issues/701)
- `oncoplot` can now summarize `variant_classifications` similar to cBioPortal style. Issue: [686](https://github.com/PoisonAlien/maftools/issues/686) 
- Added pathway support for `mafCompare()` or `clinicalEnrichment()`. Issue: [681](https://github.com/PoisonAlien/maftools/issues/681) 
- Added default title for side and topbar plots to `oncoplot`. Issue: [682](https://github.com/PoisonAlien/maftools/issues/682) 
- Added `annotationOrder` argument to `coOncoplot`. Issue: [676](https://github.com/PoisonAlien/maftools/issues/676) 
- Added `plot` argument to `survGroup`. Thank you [OmarElAshkar](https://github.com/OmarElAshkar) PR: [674](https://github.com/PoisonAlien/maftools/issues/674) 
- Added `rmFlags` argument to `read.maf`. Issue: [668](https://github.com/PoisonAlien/maftools/issues/668)
- Added `path_order` argument to `oncoplot` for custom ordering of pathways on oncoplot. 
- Added `geneMar` argument to `coBarplot`. Issue: [260](https://github.com/PoisonAlien/maftools/issues/260) 

## NEW FUNCTIONS
- `cancerhotspots` Genotype known cancer hotspots from the tumor BAM file
- `bamreadcounts` extract nucleotide counts for targeted variants from the BAM file.
- `maftools` now natively loads TCGA cohorts. `tcgaAvailable` and `tcgaLoad` will display and load the desired cohorts.
- Added `MAF` constructor function
- Added `maf2mae` for converting `MAF` to `MultiAssayExperiment` class objects Issue: [640](https://github.com/PoisonAlien/maftools/issues/640) [293](https://github.com/waldronlab/MultiAssayExperiment/pull/293) Discussion: [285](https://github.com/waldronlab/MultiAssayExperiment/discussions/285)
- Added `plotProtein` and `mafbarplot`

# CHANGES IN VERSION 2.6.05 (BioC 3.12 version)
## BUG REPORTS
- `mafSurvival` Remove non functional fn argument. Issue: [660](https://github.com/PoisonAlien/maftools/issues/660)
- `OncogenicPathways` zero entries bug and correct documentation. Issue: [656](https://github.com/PoisonAlien/maftools/issues/656)
- `plotApobecDiff` y-axis limits. Issue: [642](https://github.com/PoisonAlien/maftools/issues/642) [629](https://github.com/PoisonAlien/maftools/issues/642)
- `Rainfallplot` arrowhead bug. Issue: [628](https://github.com/PoisonAlien/maftools/issues/628) [629](https://github.com/PoisonAlien/maftools/issues/629)

## ENHANCEMENTS
- Improve `clinicalEnrichment` odd-ratio interpretation. Issue: [633](https://github.com/PoisonAlien/maftools/issues/633)
In earlier versions, odds-ratio indicated the odds of observing WT in the group of interest compared to Mutant. From this update, odds-ratio indicate the odds of observing mutant in the group of interest compared to wild-type. This way is much intuitive and easier interpret. See issue [633](https://github.com/PoisonAlien/maftools/issues/633) for details. P-values and other details are unaffected.
- Changes in how `Multi_Hit` are reported. By default two distinct types of mutations in the same gene or same patient were classified as Multi_Hit (e.g; Missense + Splice_Site = Multi_Hit; Missense + Missense = Missense). This update onward, regardless of type of mutations, if there are >1 mutations in the same gene of same patient, they are classified as Multi_Hit (e.g; Missense + Splice_Site = Multi_Hit; Missense + Missense = Multi_Hit). Issue: [347](https://github.com/PoisonAlien/maftools/issues/347) [347](https://github.com/PoisonAlien/maftools/issues/347)

# CHANGES IN VERSION 2.6.00
## BUG REPORTS
- Incorrect deduplication cases while validating MAFs. Issue: [623](https://github.com/PoisonAlien/maftools/issues/623)
- Fix repeated domain labels in `lollipopPlot2`. Issue: [614](https://github.com/PoisonAlien/maftools/issues/614)
- Fix Copy number labels and coloring in `coBarPlot`. Issue: [609](https://github.com/PoisonAlien/maftools/issues/609)
- Fix coloring issues for annotations in `oncoplot` in R < 4.0. Issue: [599](https://github.com/PoisonAlien/maftools/issues/599)

## ENHANCEMENTS
- Updated TCGA TMB. All variants from MC3 results are restricted and harmonized to Agilent 35.8 MB capture kit. See `?tcgaCompare` for details. Issue: [612](https://github.com/PoisonAlien/maftools/issues/612)
- Added `SIMC1` protein to domain database. Issue: [616](https://github.com/PoisonAlien/maftools/issues/616)
- Added `sampleOrder1` and `sampleOrder2` arguments to `coOncoplot`. Issue: [592](https://github.com/PoisonAlien/maftools/issues/592)
- Added `gene_mar` `outer_mar` argument to `coOncoplot`. Issue: [260](https://github.com/PoisonAlien/maftools/issues/260)
- Added sample size warning messages to `mafCompare`. Issue: [602](https://github.com/PoisonAlien/maftools/issues/602)
- Added `cohortFontSize` and `axisFontSize` to `tcgaCompare`

## NEW FUNCTIONS
- Added `filterMaf` function. 

## DEPRECATED
- `signatureEnrichment` will be deprecated in future. Currently kept for legacy purpose with a warning message. Issue:  [607](https://github.com/PoisonAlien/maftools/issues/607) [615](https://github.com/PoisonAlien/maftools/issues/615)

