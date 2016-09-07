# Zhian's simuPOP simulation scripts

## Purpose

The purpose of these scripts is to simulate populations with varying levels of
sexual reproduction and eventually varying levels of different population
genetic strategies such as admixture. 

## Current software environment:

 - Python 3.4
 - simuPOP 1.1.7
 - numpy 1.11.0

I used conda to install all of these. [Here's the intro to conda][conda]. It's
worth noting that I did have some trouble getting everything to work correctly
since simuPOP can't be installed with python 3.5 (for some weird reason, black
magic, maybe? I dunno). These are my steps:

1. Install conda (on my OSX, I used `brew install Caskroom/cask/anaconda`, but
   there are instructions on the site.)
2. Update my modules via conda, Downgrade python and install simuPOP:

```
conda update --all 
conda install python=3.4 
conda install -c https://conda.binstar.org/bpeng simuPOP
``` 

Finally, I made sure that numpy and feather-format were installed:

```
conda install numpy
pip install feather-format # for transfer between python and R
```

## R package

To facilitate analysis of the simulations, I've written an R package called
"zksimanalysis". It can be installed from the root directory in R with:

```r
devtools::install("zksimanalysis")
```


## Current workflow/implementation

(2016-08-18)

Currently, the script can be executed from the top-level directory. I still need
to set this up so it dumps the results in a "results" directory, but that will
come later as I need it. All of the simulation scripts will live in the 
[`simulations/`][simulations] directory.

With default values:

```sh
python simulations/simulate_sex_rate.py
```

With a configuration file called `CONFIG.args`:

```sh
python simulations/simulate_sex_rate.py @CONFIG.args
```

Once all of these are created, they can be gathered up into the feather format:

```sh
python organizing/dir2feather.py --regex gen_10 --prefix foo --group_by sex --zip --out pillow
```

## [Modules][modules]

The `modules/` directory contains python scripts that can be loaded
as modules for classes and functions.

#### [random_alleles.py][random_alleles]

This module implements two classes, 

 - *zk_locus* will initialize a locus with a mutation rate, allele frequencies,
   repeat length, and allele names.
 - *zk_loci* contains a whole bunch of *zk_locus* classes. This will provide an
   accessor to all of the loci names, allele names, allele frequencies, etc.

#### [zk_utils.py][zk_utils]

This contains misc functions for evolution


#### [pop2df.py][pop2df]

This contains the code for converting population objects to pandas data frames
with the function `pops2df()`.


## [Simulation Scripts][simulations]

The following scripts are used for the simulations. Again, they are to be run 
from the top-level directory.

#### simulate_sex_rate.py

This will simulate populations with varying levels of sexual reproduction.
Currently, there is no population structure for the populations, but this can
change.


## [Organization][organizing]

This directory contains scripts used for organizing and transferring the output.

#### [dir2feather.py][dir2feather]

This will trawl a directory and convert the populations to the [feather
format][feather], which allows crosstalk between python and R for faster
read/write than traditional csv.

[conda]: http://conda.pydata.org/docs/intro.html
[modules]: ./modules
[simulations]: ./simulations
[random_alleles]: ./modules/random_alleles.py
[zk_utils]: ./modules/zk_utils.py
[pop2df]: ./modules/pop2df.py
[dir2feather]: ./organizing/dir2feather.py
[organizing]: ./organizing
[feather]: https://blog.rstudio.org/2016/03/29/feather/