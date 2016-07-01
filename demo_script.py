#!/usr/bin/env python3.4

import simuPOP as sim
import os, sys, datetime, re, multiprocessing, random, numpy
from simuOpt import setOptions
from simuPOP import utils
from simuPOP.utils import export
from simuPOP.utils import saveCSV
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
sexytime = 0.01
nloc = 10
nall = 10
murate = 1e-5
STEPS = 1000
GENERATIONS = 500
POPSIZE = 500
sexytime = sexytime*POPSIZE
SAVEPOPS = False
# This variable is used to set up information fields that are collected at
# mating. They include:
#  - clone_proj: the average number of clonal events for each parent
#  - sex_proj: The average number of sexual events for each parent
#  - mother_id: ID of the mother (0 if no mother) [SIMUPOP PARAM]
#  - father_id: ID of father (0 if no father)     [SIMUPOP PARAM]
#  - tmsrsr: time to most recent sexual reproduction
infos = ['clone_proj', 'sex_proj', 'ind_id', 'mother_id', 'father_id', 'tsmrsr']
# Initializing a population of 100 individuals with two loci each on separate
# chromosomes. These loci each have nall alleles.
allele_names = get_allele_names(nloc, nall + 1)
loci_names = get_loci_names(nloc, allele_names)
murate = scale_mutation_rate(murate, allele_names)

pop = sim.Population(
    size = POPSIZE, 
    loci = [1]*nloc, 
    lociNames = loci_names, 
    alleleNames = allele_names,
    infoFields = infos
    )

#------------------------------------------------------------------------------#
# Generating initializing operators.
#------------------------------------------------------------------------------#

# This will generate and plot allele probabilities for each locus.
loclist = get_allele_probabilities(nloc, nall)
plot_allele_probabilities(loclist, nall, loci_names)

# This will put all of the initializing steps into a list. 
# 1. initialize sex for the populations
# 2. initialize genotypes for each locus separately.
inits = list()
inits.append(sim.InitSex())
inits.append(sim.InitInfo(0, infoFields = infos))
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
males = r"Males: %d"
het = r"Het: " + r"%.2f "*nloc
generations = r"Gen: %d"


'''
Account for the genealogical history of each isolate. This is a tagger that
updates information fields for each individual produced during mating. This is
used in during mating operations. 

Parameters: 
    clone_proj: a list containing the average number of clonal events for each
        parent.
    sex_proj: a list containing the average number of sexual events for each
        parent.
    tsmrsr: Time Since Most Recent Sexual Reproduction. This is only used if
        a clonal event occurs.
'''
def update_sex_proj(clone_proj, sex_proj, tsmrsr):
    # Sexual reproduction: average clone and sex. Add one to sex.
    if len(clone_proj) > 1:
        out_sex_proj = 1 + ((sex_proj[0] + sex_proj[1]) / 2)
        out_clone_proj = ((clone_proj[0] + clone_proj[1]) / 2)
        out_tsmrsr = 0

    # Clonal reproduction: add one to clone.
    else:
        out_sex_proj = sex_proj[0]
        out_clone_proj = clone_proj[0] + 1
        out_tsmrsr = tsmrsr[0] + 1

    return out_clone_proj, out_sex_proj, out_tsmrsr

# Joining the statistics together with pipes.
stats = " | ".join([head,popsize, males, het, generations, foot])

# Heterozygosity must be evaluate for each locus. This is a quick and dirty
# method of acheiving display of heterozygosity at each locus. 
locrange = map(str, range(nloc))
lochet = '], heteroFreq['.join(locrange)

# The string for the evaluation of the stats.
stateval = " % (popSize, numOfMales, heteroFreq["+lochet+"], gen)"

# Stat and PyEval are both classes, so they can be put into variables. These
# will be evaluated as the program runs. 
statargs = sim.Stat(popSize=True, numOfMales=True, heteroFreq = range(nloc), step = STEPS)
evalargs = sim.PyEval(stats + stateval, step = STEPS)


#------------------------------------------------------------------------------#
# Generating mating.
#------------------------------------------------------------------------------#

# There are options here are things to do during mating:
# 1. Tag the parents (which fills mother_idx and father_idx)
# 2. Count the reproductive events
mate_ops = [sim.PedigreeTagger(), sim.PyTagger(update_sex_proj)]

rand_mate = sim.RandomMating(
    numOffspring=(sim.UNIFORM_DISTRIBUTION, 1, 3),
    subPops = 0, 
    weight = sexytime, 
    ops = [
        sim.MendelianGenoTransmitter(),
        sim.PedigreeTagger(infoFields=['mother_id', 'father_id']),
        sim.PyTagger(update_sex_proj),
        sim.IdTagger()
        ]
    )

# From Bo Peng:
# When there is only one parent, PedigreeTagger only uses the first field 
# (mother_id) regardless of the sex of parent. 
# 
# Because of this, I need to reassign the sex of the parents.
def reassign_parents(pop):
    for i in pop.individuals():
        if i.sex() == 1 and i.tsmrsr > 0:
            i.father_id = i.mother_id
            i.mother_id = 0.0
    return True

clone_mate = sim.HomoMating(
    chooser = sim.RandomParentChooser(),
    generator = sim.OffspringGenerator(
        ops = [
            sim.CloneGenoTransmitter(),
            sim.PedigreeTagger(infoFields=['mother_id', 'father_id']),
            sim.PyTagger(update_sex_proj)
            ],
        numOffspring=(sim.UNIFORM_DISTRIBUTION, 1, 3)
    ),
    subPops = 0, 
    weight = POPSIZE - sexytime
    )
mate_scheme = sim.HeteroMating([rand_mate, clone_mate])




postlist = list()
postlist.append(statargs)
postlist.append(evalargs)
postlist.append(sim.PyOperator(func = reassign_parents, step = STEPS))

if SAVEPOPS is True:
    outfile = "!'gen_%d.pop' % (gen)"
    finals = sim.SavePopulation(output = outfile, step = STEPS)
    postlist.append(finals)


#------------------------------------------------------------------------------#
# Do the evolution.
#------------------------------------------------------------------------------#
sim.IdTagger().reset(1) # IdTagger must be reset before evolving.
pop.evolve(
    initOps = inits,
    matingScheme = mate_scheme,
    preOps = [
        sim.StepwiseMutator(rates = murate, loci = range(nloc)),
        sim.IdTagger(),
        ],
    postOps = postlist,
    gen = GENERATIONS
    )




moms = pop.indInfo('mother_id')
dads = pop.indInfo('father_id')
tsmrsr = pop.indInfo('tsmrsr')
sim.dump(pop)
sample = sim.sampling.drawRandomSample(pop, sizes = 100)
export(pop = sample,
        format = 'csv',
        output = "deleteme.csv", 
        infoFields = ['tsmrsr', 'ind_id', 'mother_id', 'father_id'],
        sexFormatter = {1:"M", 2:"F"}, 
        affectionFormatter = None,
        gui = False
        )



def whos_got_the_keys(gus, key):
    if key in gus:
        gus[key] += 1
    else:
        gus[key] = 1
    return gus

daddict = dict()
momdict = dict()
sexdict = dict()

# Individuals carry information of who their parents were. If they were clonally
# produced, then one of the parental indices will have a -1, indicating that
# there was no contribution from that parental unit. This loops through and
# assesses what clones came from mothers and what clones came from fathers.
for i in range(POPSIZE):
    sexdict = whos_got_the_keys(sexdict, tsmrsr[i])
    if moms[i] <= 0:
        #print("dad " + str(dads[i]))
        daddict = whos_got_the_keys(daddict, dads[i])
    elif dads[i] <= 0:
        #print("mom " + str(moms[i]))
        momdict = whos_got_the_keys(momdict, moms[i])


def print_vals(indict, vals, counts):

    padding = abs(len(vals) - len(counts))
    valdict = dict()
    for i in indict.keys():
        valdict = whos_got_the_keys(valdict, indict[i])

    sys.stdout.write(vals + ":" + " "*padding + "\t")
    for i in valdict.keys():
        sys.stdout.write(str(i) + "\t")
    print("")

    sys.stdout.write(counts + ":" + " "*padding + "\t")
    for i in valdict.keys():
        sys.stdout.write(str(valdict[i]) + "\t")
    print("")

print("| "+"="*80)
print("| ______ Dads: " + str(len(daddict)))
print_vals(daddict, "| Clones Produced", "| Clonefathers")
print("| "+"-"*80)
print("| ______ Moms: " + str(len(momdict)))
print_vals(momdict, "| Clones Produced", "| Clonemothers")
print("| "+"="*80)
print("Amount of Sex: " + str(sexytime) + " individuals per generation.")
print("Time Since Sex:")
print("\tTime\tCount")
for i in sorted(sexdict.keys()):
    print("\t" + str(i) + "\t" + str(sexdict[i]))

