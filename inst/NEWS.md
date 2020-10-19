# CHANGES IN VERSION 2.6.00 (For upcoming BioC 3.12 release)
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

