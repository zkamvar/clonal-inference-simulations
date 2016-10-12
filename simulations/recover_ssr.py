#!/usr/bin/env python3

# Useful for importing any python modules I have in the "modules" directory.
# http://stackoverflow.com/a/35259170/2752888
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'modules'))
import random_alleles as ra
from zk_utils import *
import argparse
import inspect

parser = argparse.ArgumentParser(
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    fromfile_prefix_chars = '@'
    )

parser.convert_arg_line_to_args = convert_arg_line_to_args

parser.add_argument(
    "--murate", 
    type = float,
    nargs = "+",
    help = "mutation rate per locus",
    default = [1e-05]*10
    )
parser.add_argument(
    "--popfile",
    type = str,
    help = "relative path to last saved simulation pop file."
    )
# Setting up options
import simuOpt
simuOpt.setOptions(optimized = True, 
    gui = False, 
    debug = 'DBG_WARNING',
    alleleType = 'long', 
    quiet = False, 
    numThreads = 0)
import simuPOP as sim
from simuPOP.utils import export
from simuPOP.utils import saveCSV
from simuPOP.sampling import drawRandomSample
from simuPOP.sampling import drawRandomSamples

'''
A wrapper to generate a mixed mating scheme with a rate of sexual reproduction.

Parameters:
    sexrate a rate of sexual reproduction between 0 and 1

Output:
    a simuPOP Mating scheme

Examples:
    mix_mating(0.5)
'''
def mix_mating(sexrate):
    # Mating Schemes
    # ==========================================================================
    # There are options here are things to do during mating:
    # 1. Tag the parents (which fills mother_idx and father_idx)
    # 2. Count the reproductive events

    rand_mate = sim.RandomMating(
        numOffspring=(sim.UNIFORM_DISTRIBUTION, 1, 3),
        ops = [
            # sim.Stat(numOfMales = True, popSize = True),
            # sim.IfElse('numOfMales == popSize or numOfMales == 0',
            #     sim.CloneGenoTransmitter(),
            #     sim.MendelianGenoTransmitter()
            #     ),
            sim.MendelianGenoTransmitter(),
            sim.PedigreeTagger(infoFields=['mother_id', 'father_id']),
            sim.PyTagger(update_sex_proj),
            sim.IdTagger()
            ],
        subPops = 0, 
        weight = sexrate
        )

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
        weight = 1 - sexrate
        )
    return(sim.HeteroMating([rand_mate, clone_mate]))



# ==============================================================================
# MAIN PROGRAM
# ==============================================================================

if __name__ == '__main__':
    pars   = parser.parse_args()
    infile = os.path.abspath(pars.popfile)
    the_path = os.path.dirname(infile)
    the_dataset = os.path.basename(infile)
    a, seed, b, sex, c, gen, d, rep = the_dataset.split(".p")[0].split("_")
    # seed_0_sex_0.0001_gen_10000_rep_06.pop
    GENRATIONS = 10000 - int(gen)
    sexrate = float(sex)
    cwd = os.getcwd()
    os.chdir(the_path)
    pop  = sim.loadPopulation(the_dataset)
    sim.IdTagger().reset(1) # IdTagger must be reset before evolving.
    nloc = sum(list(pop.numLoci()))
    clonemate = sim.HomoMating(
        chooser = sim.RandomParentChooser(),
        generator = sim.OffspringGenerator(
            ops = [
                sim.CloneGenoTransmitter(),
                sim.PedigreeTagger(infoFields=['mother_id', 'father_id']),
                sim.PyTagger(update_sex_proj)
                ],
            numOffspring = (sim.UNIFORM_DISTRIBUTION, 1, 3)
            )
        )
    EVOL = dict()
    STEPS = 1000
    head, foot  = r"'", r"\n'"
    popsize     = r"Pop Size: {}"
    males       = r"Males: {}"
    generations = r"Gen: {:05d}"
    # Joining the statistics together with pipes.
    stats = " | ".join([head, popsize, males, generations, foot])

    # Heterozygosity must be evaluate for each locus. This is a quick and dirty
    # method of acheiving display of heterozygosity at each locus. 
    # locrange = map(str, range(nloc))
    # lochet = '], heteroFreq['.join(locrange)


    # The string for the evaluation of the stats.
    # stateval = " % (popSize, numOfMales, heteroFreq["+lochet+"], gen, rep)"
    # stateval  = " % "
    stateval = ".format("
    stateval += "popSize, "
    stateval += "numOfMales, "
    stateval += "gen"
    stateval += ")"

    # Stat and PyEval are both classes, so they can be put into variables. These
    # will be evaluated as the program runs. 
    statargs = sim.Stat(
        popSize = True, 
        numOfMales = True,
        step = STEPS
    )
    evalargs = sim.PyEval(stats + stateval, step = STEPS)

    EVOL['postOps'] = [
        statargs,
        evalargs#,
        # sim.PyOperator(func = reassign_parents, step = STEPS),
    ]
    seedf             = "seed_" + seed
    sexf              = "_sex_" + sex
    sexseed           = seedf + sexf
    outfile           = "!'"+ sexseed + "_gen_{:05d}_rep_"+rep+".pop'"
    outfile           = outfile + ".format(gen)"
    EVOL['postOps']  += [sim.SavePopulation(output = outfile, step = STEPS)]
    
    GENERATIONS = 10000 - int(gen)
    print(GENERATIONS)
    EVOL['gen'] = GENERATIONS
    preclone = [
        sim.Stat(numOfMales = True, popSize = True),
        sim.StepwiseMutator(rates = pars.murate, loci = range(nloc)),
        sim.IdTagger(),
    ]
    
    if sexrate == 0.0:
        EVOL['matingScheme'] = clonemate
        EVOL['preOps'] = preclone
    else:
        EVOL['matingScheme'] = mix_mating(sexrate)
        EVOL['preOps'] = [
            sim.Stat(numOfMales = True, popSize = True),
            sim.TerminateIf('numOfMales == 0 or numOfMales == popSize'),
            sim.StepwiseMutator(rates = pars.murate, loci = range(nloc)),
            sim.IdTagger(),
        ]
    res = pop.evolve(**EVOL)
    if (res < GENERATIONS):
        EVOL['matingScheme'] = clonemate
        EVOL['preOps'] = preclone
        EVOL['gen'] = GENERATIONS - res
        pop.evolve(**EVOL)
    os.chdir(cwd)
