# CHANGES IN VERSION 2.10.01

- Added PALB2 protein structure. Issue: [677](https://github.com/PoisonAlien/maftools/issues/677)

## BUG FIXES
- `annovarToMaf` struggles with Variant_Type. Issue: [803](https://github.com/PoisonAlien/maftools/issues/803)
- `lollipopPlot2` and Domain Details for Further Analysis Issue: [794](https://github.com/PoisonAlien/maftools/issues/794)
- `coOncoplot` Ignores Annotation Colours. Issue: [786](https://github.com/PoisonAlien/maftools/issues/786)
- `subsetMaf` Lost Columns if clinQuery Result is One Row. Issue: [785](https://github.com/PoisonAlien/maftools/issues/785)

# CHANGES IN VERSION 2.10.00 (BC 3.14 release version)

## NEW FUNCTIONS
- `sampleSwaps` Given a list BAM files, the function genotypes known SNPs and identifies potentially related samples.

## BUG FIXES
- Return `mutCountMatrix` output as a matrix Issue: [769](https://github.com/PoisonAlien/maftools/issues/769)

## ENHANCEMENTS
- Added `showPct` argument to `oncoplot`. Issue: Issue: [771](https://github.com/PoisonAlien/maftools/issues/780)
- Silently return sample order from `oncoplot` and `coOncoplots` Issue: [771](https://github.com/PoisonAlien/maftools/issues/771)
