# zksimanalysis 0.10.0.9000

## Dependency updates

* Non-cran repositories have been removed and minimum versions for poppr and 
  adegenet have been updated

# zksimanalysis 0.9.0.9000

## NEW FUNCTIONS

* `tidy_jackia()` will perform analysis of `ia()` on data reduced to the number
  of observed MLGs. This will give an estimate of the spread of the data at that
  sample size.

# zksimanalysis 0.8.1.9000

## NEW FEATURES

* `roc()` gains the argument `flip` which specifies a column of logicals to use
  that would flip a negative decision positive. 

# zksimanalysis 0.8.0.9000

## NEW FUNCTION

* `filter_maf()` will filter genlight objects based on minor allele frequency.

## NEW FEATURES

* `shuffle_genlight()` can now shuffle blocks of SNPs.

# zksimanalysis 0.7.0.9000

## NEW FUNCTION

* `roc()` calculates the receiver operator characteristic.

## NEW DEPENDENCIES

* The package _tidyverse_ has been added to replace _dplyr_ and _magrittr_. This
  also adds _purrr_.
* The color package _viridis_ has been added for plotting.
* The _flux_ package has been added for calculation of area under ROC curves.

# zksimanalysis 0.6.0.9000

## NEW FUNCTION

* The manipulation functions `pop_NA()` and `pop_mutator()` have been exported.

## NEW FEATURES

* `tidy_ia()` can now take genomic data
* `genomic_ia()` now has a reporter when interactive.
* parallel processing is hard-coded to FALSE in shuffling and data import for
  genomic data.

# zksimanalysis 0.5.0.9000

## NEW FUNCTION

* `feather2genlight()` will convert feather formatted files to genlight/snpclone
  formatted data. 
* `genomic_ia()` Calculates and resamples the standardized index of association
  from your data.
* `sample_bitwise_ia()` Resamples your data once and calculates the standardized
  index of association
* `shuffle_genlight()` shuffles a genlight/snpclone object to unlink markers.

## MISC

* Added internal functions from
  https://github.com/grunwaldlab/supplementary-poppr-2.0/blob/master/Rscripts/my_functions.R

# zksimanalysis 0.4.0.9000

## NEW FUNCTION

* `tidy_locus_contribution()` will calculate pgen, genotype accumulation curve,
  round-robin multilocus genotypes, and the diversity of contributions of each
  locus to the number of multilocus genotypes.

# zksimanalysis 0.3.1.9000

## ENHANCEMENTS

* The function `tidy_locus_table()` will now report clone-corrected data.
* Clone-corrected data from `tidy_ia()` will now undergo permutation testing.

# zksimanalysis 0.3.0.9000

## NEW FUNCTION

* `tidy_locus_table()` is here to calculate the diversity of your loci.

# zksimanalysis 0.2.0.9000

## NEW FUNCTION

* The new function `tidy_diversity()` will get diversity estimates and return
  a tidy data frame with the estimates and their variances.
  
## NEW IMPORTS

* poweRlaw is now imported.

# zksimanalysis 0.1.1.9000

* The `sample` argument for `feather2genind()` can now be NULL to read in the
  entire file.
* The function `tidy_ia()` will now clone-correct data.

# zksimanalysis 0.1.0.9000

* Added a `NEWS.md` file to track changes to the package.
* add new function `listfun()` to process single-item list colums that can be
  stored in dplyr.

## BUG FIX

* A bug with `feather2genind()` where `pop` was being added in as a locus was
  fixed.

# zksimanalysis 0.0.0.9000

* Initial package version
