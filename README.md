# Zhian's simuPOP simulation scripts

## Purpose

The purpose of these scripts is to simulate populations with varying levels of
sexual reproduction and eventually varying levels of different popualtion
genetic strategies such as admixture. 

## Current software environment:

 - Python 3.4
 - simuPOP 1.1.7
 - numpy 1.11.0

I used conda to install all of these. [Here's the intro to
conda](http://conda.pydata.org/docs/intro.html). It's worth noting that I did
have some trouble getting everything to work correctly since simuPOP can't be
installed with python 3.5 (for some weird reason, black magic, maybe? I dunno).
These are my steps:

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

## Original workflow/implementation

The original way this was done was to generate configuration files with
`CFGenerator.py` and then to use `Rand_Batches.pl` to first generate all of the
burnin populations. The point of the burnins was to have populations undergo
sexual recombination for 1000 generations to ensure equilibrium. After the
burnins were generated, all of the simulations would begin under the specified
parameters.

## File Descriptions

- `random_alleles.py` contains a script to generate random allele frequencies for
a given number of loci. This is based off of normalized draws from a uniform
distribution. 
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
- `demo_script.py` is a script that is in progress that hard codes all of the
parameters in an attempt to understand what is going on. 