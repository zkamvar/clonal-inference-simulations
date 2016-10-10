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
