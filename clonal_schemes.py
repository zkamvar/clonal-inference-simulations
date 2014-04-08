#!/usr/bin/env python

import simuPOP as sim
import os, sys, datetime, commands, re, exceptions, multiprocessing, random, numpy
from simuOpt import setOptions
from simuPOP import utils
from simuPOP.utils import export
from simuPOP.sampling import drawRandomSample, drawRandomSamples 
from random_alleles import *


#------------------------------------------------------------------------------#
# clone simulator in case the population goes strictly clonal.
# newG0 is the number of generations left for the population to go.
# This will also create a population marking the point at which all individuals
# are of the same sex (or mating type). 
#------------------------------------------------------------------------------#
def goToClone(pops, L0, newG0, G0, N0, R0):
    simu = sim.Simulator(pops)
    print("GOING TO CLONAL REPRODUCTION")
    popsize = r"'Pop Size: %d'"
    males = r"'Number of Males: %d'"
    het = r"'Heterozygosity: %.2f'"
    generations = r"'Gen: %d'"
    stats = " | ".join([popsize, males, het, generations])
    simu.evolve(  
        matingScheme = sim.RandomSelection(subPops=0),
        preOps=[sim.StepwiseMutator(rates=1e-5, loci=L0)],
        postOps=[
            # This outputs more informative statistics to the screen. 
            sim.Stat(popSize=True, numOfMales=True, heteroFreq=L0, step=G0/10),
            sim.PyEval(r"'Pop Size: %d | Number of Males: %d | Heterozygosity: \
%.2f | Gen: %d\n' % (popSize, numOfMales, heteroFreq[1], gen-("+str(G0)+"/10))", step=G0/10),
            sim.SavePopulation(output="!'"+str(N0)+"/"+str(N0)+
            "_gen_%d.pop' % (gen-("+str(G0)+"/10))", step=G0/10),
        ],
        finalOps=[
            sim.SavePopulation(output="!'"+str(N0)+"/"+str(N0)+
            "_gen_%d.pop' % (gen-("+str(G0)+"/10))", step=G0/10),
        ],
        gen = newG0
    )
    print("\n---\nRep "+str(R0)+" is done!\n---\n")


print("Hey there!")

#------------------------------------------------------------------------------#
# Variables to set up 
#------------------------------------------------------------------------------#
C0 = 90
STEPS = 10
GENERATIONS = 1000

# Initializing a population of 100 individuals with two loci each on separate
# chromosomes. These loci each have 4 alleles.
pop = sim.Population(
    size = 100, 
    loci = [1, 1], 
    lociNames = 'L1 L2'.split(), 
    alleleNames = get_allele_names(2, 4)
    )

#------------------------------------------------------------------------------#
# Generating initializing operators.
#------------------------------------------------------------------------------#

# This will generate and plot allele probabilities for each locus.
loclist = get_allele_probabilities(2, 4)
plot_allele_probabilities(loclist, 4)

# This will put all of the initializing steps into a list. 
# 1. initialize sex for the populations
# 2. initialize genotypes for each locus separately.
inits = list()
inits.append(sim.InitSex())
for i in range(len(loclist)):
    inits.append(sim.InitGenotype(freq = loclist[i], loci = i))

#------------------------------------------------------------------------------#
# Generating stats operators.
#------------------------------------------------------------------------------#
# The stats operators must be raw strings. Since it would be advantageous later
# to be able to add in different stats, it's better to join a bunch of small
# strings than write one long one.

# The stats must have single quotes around them.
foot = r"\n'"
head = r"'"

popsize = r"Pop Size: %d"
males = r"Number of Males: %d"
het = r"Heterozygosity: %.2f %.2f"
generations = r"Gen: %d"

# Joining the statistics together with pipes.
stats = " | ".join([head,popsize, males, het, generations, foot])

# The string for the evaluation of the stats. 
stateval = " % (popSize, numOfMales, heteroFreq[0], heteroFreq[1], gen)"

# Stat and PyEval are both classes, so they can be put into variables. These
# will be evaluated as the program runs. 
statargs = sim.Stat(popSize=True, numOfMales=True, heteroFreq=[0, 1], step = STEPS)
evalargs = sim.PyEval(stats + stateval, step = STEPS)


#------------------------------------------------------------------------------#
# Generating mating.
#------------------------------------------------------------------------------#

rand_mate = sim.RandomMating(subPops = 0, weight = C0)
clone_mate = sim.RandomSelection(subPops = 0, weight = 100 - C0)
mate_scheme = sim.HeteroMating([rand_mate, clone_mate])
outfile = "!'gen_%d.pop' % (gen)"
finals = sim.SavePopulation(output = outfile, step = STEPS)
#------------------------------------------------------------------------------#
# Do the evolution.
#------------------------------------------------------------------------------#
pop.evolve(
    initOps = inits,
    matingScheme = mate_scheme,
    preOps = [sim.StepwiseMutator(rates = 1e-5, loci = [0,1])],
    postOps = [statargs, evalargs, finals],
    gen = GENERATIONS
    )
sim.dump(pop)



