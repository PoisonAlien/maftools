# CHANGES IN VERSION 2.21.2

## BUG FIXES
- Bug fix in using `keepGeneOrder` in `coOncoplot()`. Issue: [1061](https://github.com/PoisonAlien/maftools/issues/1061)
- Bug fix in using `selectedPathways` in `oncoplot()`. Issue: [1041](https://github.com/PoisonAlien/maftools/issues/1041)
- Add an error message when bai files are missing `sampleSwaps()`. Issue: [1028](https://github.com/PoisonAlien/maftools/issues/1028)
- Bug fix in `tmb` while handling multiple MAFs. Issue: [1018](https://github.com/PoisonAlien/maftools/issues/1018)
- Handle missing `NA`s while sub-setting for ranges. Issue: [1013](https://github.com/PoisonAlien/maftools/issues/1013)
- Better error handling when zero mutated samples are encountered in `clinicalEnrichment`. Issue: [1010](https://github.com/PoisonAlien/maftools/issues/1010)
- MAJOR: `read.maf` by default coerces clinical data columns to character. This bug fix avoids it and is auto detected. Issue: [997](https://github.com/PoisonAlien/maftools/issues/997)

## ENHANCEMENTS
- Better handling of color codes for continuous variable annotations [1053](https://github.com/PoisonAlien/maftools/issues/1053)
- Add `right_mar` to `gisticOncoPlot`[1043](https://github.com/PoisonAlien/maftools/issues/1043)
- Added `PPDPFL` to protein domain database Issue: [1025](https://github.com/PoisonAlien/maftools/issues/1025)
- Better sorting of oncoplot with `collapsePathway`
- Changed default background for oncoplot from `gray` to `#ecf0f1`
- Changed default signature database to SBS_v3.4 (from legacy)
- Update `tmb` function

## BREAKING CHNAGES
- Column order required for `pathways()` function changed from Pathway,Gene to Gene,Pathway. Issue: [1041](https://github.com/PoisonAlien/maftools/issues/1041)

## NEW FUNCTIONS
- `gisticCompare()` for comparing two GISTIC objects
- `segSummarize()` for summarizing DNAcopy segments
