# Zhian's simuPOP simulation scripts

## Purpose

The purpose of these scripts is to simulate populations with varying levels of
sexual reproduction and eventually varying levels of different popualtion
genetic strategies such as admixture. They rely on simuPOP and numpy for python
v. 2.7. All of these files were modified last in early 2013 and need to undergo
heavy revamping to make them portable.

## Original workflow/implementation

The original way this was done was to generate configuration files with `CFGenerator.py`
and then to use `Rand_Batches.pl` to first generate all of the burnin populations.
The point of the burnins was to have populations undergo sexual recombination for
1000 generations to ensure equilibrium. After the burnins were generated, all of the
simulations would begin under the specified parameters.

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