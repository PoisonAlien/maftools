# CHANGES IN VERSION 3.0.00 (BioC 3.14 RC)

## NEW FUNCTIONS
- `sampleSwaps` Given a list BAM files, the function genotypes known SNPs and identifies potentially related samples.

## BUG FIXES
- Return `mutCountMatrix` output as a matrix Issue: [769](https://github.com/PoisonAlien/maftools/issues/769)

## ENHANCEMENTS
- Added `showPct` argument to `oncoplot`. Issue: Issue: [771](https://github.com/PoisonAlien/maftools/issues/780)
- Silently return sample order from `oncoplot` and `coOncoplots` Issue: [771](https://github.com/PoisonAlien/maftools/issues/771)
