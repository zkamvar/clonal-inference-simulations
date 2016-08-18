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

Finally, I made sure that numpy was installed:

```
conda install numpy
```

## Current workflow/implementation

(2016-08-18)

Currently, the script can be executed from the top-level directory. I still need
to set this up so it dumps the results in a "results" directory, but that will
come later as I need it. All of the simulation scripts will live in the 
[`simulations/`][simulations] directory.

With no configuration file (options will be presented and the configuration saved):

```sh
python simulations/simulate_sex_rate.py
```

With a configuration file called `config.cfg`:

```sh
python simulations/simulate_sex_rate.py --config config.cfg
```

## Modules

The [`modules/`][modules] directory

#### random_alleles.py

This module implements two classes, 

 - *zk_locus* will initialize a locus with a mutation rate, allele frequencies,
   repeat length, and allele names.
 - *zk_loci* contains a whole bunch of *zk_locus* classes. This will provide an
   accessor to all of the loci names, allele names, allele frequencies, etc.

## Simulation Scripts

The following scripts are used for the simulations. Again, they are to be run 
from the top-level directory.

#### simulate_sex_rate.py

This will simulate populations with varying levels of sexual reproduction.
Currently, there is no population structure for the populations, but this can
change.

## Old File Descriptions


- `CFGenerator.py` will generate configuration files for the script. 
- `multi_burnin.py` will generate named burnin files.
- `burnin.py` will generate a burnin files named BURNIN.pop
- `Sim_Rand_Seed.py` will simulate a population from a random seed population
generated from `CFGenerator.py` and `multi_burnin.py`.
- `Sim_Single_Seed.py` will simulate a population from a seed titled `BURNIN.pop`
generated from `burnin.py`
- `Rand_Batches.pl` will perform burnins and then run all simulations in batch. 
- `gather_finished.sh` will gather all of the finished simulation files into a
single directory.
- `demo_script.py` serves as scratch space


## Original workflow/implementation

The original way this was done was to generate configuration files with
`CFGenerator.py` and then to use `Rand_Batches.pl` to first generate all of the
burnin populations. The point of the burnins was to have populations undergo
sexual recombination for 1000 generations to ensure equilibrium. After the
burnins were generated, all of the simulations would begin under the specified
parameters.


[conda]: http://conda.pydata.org/docs/intro.html
[modules]: ./modules
[simulations]: ./simulations