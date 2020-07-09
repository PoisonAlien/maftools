# CHANGES IN VERSION 2.4.10 (In development)
## BUG FIXES
- Clinical features with same labels are colored similarly in Oncoplot. Issue: [567](https://github.com/PoisonAlien/maftools/issues/567)
- Legend Ordering Not Intuitive. Issue: [424](https://github.com/PoisonAlien/maftools/issues/424)

## ENHANCEMENTS
- Added `minMut` argument to `oncoplot`. Issue: [549](https://github.com/PoisonAlien/maftools/issues/549)
- Added `barcodeSrt`, `annoBorderCol` arguments to `oncoplot`. Issue: [564](https://github.com/PoisonAlien/maftools/issues/564) [565](https://github.com/PoisonAlien/maftools/issues/565)

# CHANGES IN VERSION 2.4.05 (Release/Bioconductor version)
## BUG FIXES
- Bugfix in annotation legend order in `oncoplot`. Issue: [541](https://github.com/PoisonAlien/maftools/issues/541)
- Error in `plotCBSsegments`. Issue: [539](https://github.com/PoisonAlien/maftools/issues/539)
- `merge_mafs` bug fix due to discrepancy in column names. Issue: [538](https://github.com/PoisonAlien/maftools/issues/538)
- `readGistic` TCGA sample names bug fix. Issue: [536](https://github.com/PoisonAlien/maftools/issues/536)
- Highlight most significant entry when a variants overlaps with multiple CNV regions Issue: [524](https://github.com/PoisonAlien/maftools/issues/524) 
- `genesToBarcodes` sample selection fix
- `coOncoplot` bug fix with `showSampleNames` argument. Issue: [306](https://github.com/PoisonAlien/maftools/issues/306)
- `removeNonMutated` bug fix with `oncoplot` argument. Issue: [530](https://github.com/PoisonAlien/maftools/issues/530)

## ENHANCEMENTS
- Make `flags` internal function with top 100 genes. Issue: [542](https://github.com/PoisonAlien/maftools/issues/542)
- Add legend to `gisticBubblePlot`. Issue: [544](https://github.com/PoisonAlien/maftools/issues/544)
- Add `bgCol` and `borderCol` to  `gisticOncoPlot`. Issue: [#533](https://github.com/PoisonAlien/maftools/issues/533)
- Allow multiple CN types in oncoplot Issue: [#490](https://github.com/PoisonAlien/maftools/issues/490)
- Allow `getSampleSummary` to include non mutated samples Issue: [#527](https://github.com/PoisonAlien/maftools/issues/527)
- Added `color` argument to `plotVaf` for manual color codes Issue: [#531](https://github.com/PoisonAlien/maftools/issues/531)
- Added `rm_zero` argument to `tcgaCompare`

# CHANGES IN VERSION 2.4.00
## BUG FIXES

- `lollipopPlot2` error with zero mutation Issue: [#480](https://github.com/PoisonAlien/maftools/issues/480)
- `inferHeterogeneity` PR [#473](https://github.com/PoisonAlien/maftools/issues/473)
- `oncoplot` cnv crosses border [#472](https://github.com/PoisonAlien/maftools/issues/472)
- Avoid wrongly naming columns in `annovarToMaf` Issue: [#457](https://github.com/PoisonAlien/maftools/issues/457)
- Duplicated color code Issue: [#452](https://github.com/PoisonAlien/maftools/issues/452)
- oncoplot with exprsTbl leads to wrong alignment Issue: [#451](https://github.com/PoisonAlien/maftools/issues/451)
- plotMafSummary with single sample Issue: [#449](https://github.com/PoisonAlien/maftools/issues/449)
- Avoid empty spaces and quoates while reading funcotator MAF Issue: [#403](https://github.com/PoisonAlien/maftools/issues/403) [#397](https://github.com/PoisonAlien/maftools/issues/397)

## ENHANCEMENTS

- `OncogenicPathways` improvements Issue: [#505](https://github.com/PoisonAlien/maftools/pull/509)
- Added `showSum`, `colNC`, `nShiftSymbols`, `sigSymbolsSize`, `sigSymbolsFontSize`, `pvSymbols` arguments to `somaticInteractions` PR: [#505](https://github.com/PoisonAlien/maftools/pull/505) Thanks [zmiimz](https://github.com/zmiimz)
- Added `anno_height`, `drawBox`, `drawExtraBar` arguments to `oncoplot` PR: [#501](https://github.com/PoisonAlien/maftools/pull/501) Thanks [Kai Gu](https://github.com/kaigu1990)
- Added `decreasing` argument to `tcgaComapre` Issue: [#497](https://github.com/PoisonAlien/maftools/issues/497)
- Added startgain and startlost to `annovarToMaf` Issue: [#487](https://github.com/PoisonAlien/maftools/issues/487)
- `y_lims` argument in gisticChromPlot Issue: [#485](https://github.com/PoisonAlien/maftools/issues/485)
- Added support for highlighting multiple additionalFeatures Issue: [#476](https://github.com/PoisonAlien/maftools/issues/476)
- Pairwise t-test and mutational load in `tcgaCompare` Issue: [#453](https://github.com/PoisonAlien/maftools/issues/453)
- Added `titleText` argument to `oncoplot` Issue: [#448](https://github.com/PoisonAlien/maftools/issues/448)
- Highlight mutated genes on gisticChromplot Issue: [#443](https://github.com/PoisonAlien/maftools/issues/443)
- Added `legend_height` argument to `oncoplot` Issue: [#346](https://github.com/PoisonAlien/maftools/issues/346)
- Added `bgBorderCol` `domainBorderCol` `showLegend` argument to `lollipopPlot`

## NEW FUNCTIONS AND FEATURES

- `tmb` Simple tumor mutation burden estimation function.
- `setdiffMAF` and `intersectMAF`. Set operations for MAF files. PR: [#508](https://github.com/PoisonAlien/maftools/pull/508) Thanks [Shixiang Wang](https://github.com/ShixiangWang)
- `tcgaDriverBP` Compare genes to known TCGA drivers and their biological pathways
- `vafCompare` plots vaf distribution of target genes from two cohorts Issue: [#495](https://github.com/PoisonAlien/maftools/issues/495)
- `coBarplot` side-by-side barplot for maf comparison Issue: [#486](https://github.com/PoisonAlien/maftools/issues/486)
- Group genes in `oncoplot` by specific pathways. This can be done by setting `pathways = 'auto'` or by providing a two column data.frame/tsv-file with Gene names and their corresponding pathway belongings.
- Oncoplot gains `rightBarData` and `leftBarData` arguments for user specific side bar plots.
- Default flatUI colors for lollipop plots

## Deprecated
- `geneCloud` and `pancanComparison`
