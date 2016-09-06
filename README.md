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

To get the options, use the `-h` flag:

```sh
$ python simulations/simulate_sex_rate.py -h
simuPOP Version 1.1.7 : Copyright (c) 2004-2016 Bo Peng
Revision 5000 (Jan 21 2016) for Python 3.4.4 (64bit, 0thread)
Random Number Generator is set to mt19937 with random seed 0x16c0c9b8ac36874a.
This is the optimized long allele version with 18446744073709551616 maximum allelic states.
For more information, please visit http://simupop.sourceforge.net,
or email simupop-list@lists.sourceforge.net (subscription required).
Turn on debug 'DBG_WARNING'
usage: simulate_sex_rate.py [-h] [--POPSIZE POPSIZE] [--nloc NLOC]
                            [--outfile OUTFILE] [--cfg CFG]
                            [--GENERATIONS GENERATIONS] [--STEPS STEPS]
                            [--sexrate SEXRATE [SEXRATE ...]]
                            [--murate MURATE [MURATE ...]] [--amin AMIN]
                            [--amax AMAX] [--rep REP] [--nseed NSEED]

optional arguments:
  -h, --help            show this help message and exit
  --POPSIZE POPSIZE     Set the census populations size. (default: 1000)
  --nloc NLOC           Set the number of unlinked loci. (default: 10)
  --outfile OUTFILE     Set the name of the output file. (default: foo)
  --cfg CFG             Set the name of the configuration file. (default:
                        CONFIG.args)
  --GENERATIONS GENERATIONS
                        Number of generations to evolve. (default: 10001)
  --STEPS STEPS         Steps at which to save evolving populations (default:
                        1000)
  --sexrate SEXRATE [SEXRATE ...]
                        Percentage of sexual reproduction (default: [0.0,
                        0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5,
                        1.0])
  --murate MURATE [MURATE ...]
                        mutation rate per locus (default: [1e-05, 1e-05,
                        1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05,
                        1e-05])
  --amin AMIN           Minimum number of alleles per locus (default: 6)
  --amax AMAX           Maximum number of alleles per locus (default: 10)
  --rep REP             Number of replicates per population (default: 10)
  --nseed NSEED         Number of seed populations. Multiply this with rep and
                        sexrate to calculate the number of populations to be
                        evolved. (default: 100)
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

#### [pop2feather.py][pop2feather]

This will convert the populations to the [feather format][feather], which allows
crosstalk between python and R for faster read/write than traditional csv.

[conda]: http://conda.pydata.org/docs/intro.html
[modules]: ./modules
[simulations]: ./simulations
[random_alleles]: ./modules/random_alleles.py
[zk_utils]: ./modules/zk_utils.py
[pop2df]: ./modules/pop2df.py
[pop2feather]: ./organizing/pop2feather.py
[organizing]: ./organizing
[feather]: https://blog.rstudio.org/2016/03/29/feather/