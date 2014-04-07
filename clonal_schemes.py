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

pop = sim.Population(
    size = 100, 
    loci = [1, 1], 
    lociNames = 'L1 L2'.split(), 
    alleleNames = '201 202 204 206'.split()
    )

loclist = get_allele_probabilities(2, 4)
plot_allele_probabilities(loclist, 4)

inits = list()
inits.append(sim.InitSex())
for i in range(len(loclist)):
    inits.append(sim.InitGenotype(freq = loclist[i], loci = i))

popsize = r"'Pop Size: %d"
males = r"Number of Males: %d\n'"
#het = r"Heterozygosity: %.2f %.2f\n'"
#generations = r"Gen: %d"
stats = r" | ".join([popsize, males])#, het, generations])
stateval = " % (popSize, numOfMales)"#, heteroFreq[0], heteroFreq[1], gen)"
stats = stats + stateval
statargs = sim.Stat(popSize=True, numOfMales=True, heteroFreq=[1, 1], step = 1)
evalargs = sim.PyEval(stats, step = 1)

pop.evolve(
    initOps = inits,
    matingScheme = sim.RandomMating(),
    preOps = [sim.StepwiseMutator(rates = 1e-5, loci = [1,1])],
    postOps = [statargs, evalargs],
    gen = 10
    )




