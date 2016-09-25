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
    "--POPSIZE", 
    type = int,
    help = "Set the census populations size.",
    default = 1000
    )
parser.add_argument(
    "--nloc",
    type = int,
    help = "Number of loci in the genome",
    default = 10000
    )
parser.add_argument(
    "--nchrom",
    type = int,
    help = "Number of chromosomes in the genome",
    default = 10
    )
parser.add_argument(
    "--outfile", 
    type = str,
    help = "Set the name of the output file.",
    default = "foo"
    )
parser.add_argument(
    "--cfg", 
    type = str,
    help = "Set the name of the configuration file.",
    default = "CONFIG.args"
    )
parser.add_argument(
    "--GENERATIONS", 
    type = int,
    help = "Number of generations to evolve.",
    default = 10001
    )
parser.add_argument(
    "--STEPS", 
    type = int,
    help = "Steps at which to save evolving populations",
    default = 1000
    )
parser.add_argument(
    "--sexrate", 
    type = float,
    nargs = "+",
    help = "Percentage of sexual reproduction",
    default = [0.0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0]
    )
parser.add_argument(
    "--murate", 
    type = float,
    help = "mutation rate per locus",
    default = 1e-05
    )
parser.add_argument(
    "--rep", 
    type = int,
    help = "Number of replicates per population",
    default = 10
    )
parser.add_argument(
    "--sample", 
    type = int,
    nargs = "+",
    help = "Sample sizes to take from the final population",
    default = [10, 25, 50, 100]
    )
parser.add_argument(
    "--nseed", 
    type = int,
    help = "Number of seed populations. Multiply this with rep and sexrate " +
    "to calculate the number of populations to be evolved. ",
    default = 100
    )
# Setting up options
import simuOpt
simuOpt.setOptions(optimized = False, 
    gui = False, 
    debug = 'DBG_WARNING',
    alleleType = 'binary', 
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
            sim.Recombinator(rates = 0.01),
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

def saveSampleReps(the_simulator, generations, maxgen, samp, outfile, sampstring):
    nrep = range(the_simulator.numRep())
    for pop, gen, rep in zip(the_simulator.populations(), generations, nrep):
        # Ensuring that we draw without replacement
        if gen < maxgen:
            continue
        p = drawRandomSample(pop, sum(s))
        current = 0
        for s in samp:
            current += s
            outfile = outfile + "_sam_{:"+sampstring+"d}.pop"
            outfile = outfile.format(gen, rep, s)
            out = p.extractIndividuals(indexes = range(current, s))
            out.SavePopulation(outfile)


def sim_partial_clone(loci, sexrate, STEPS, GENERATIONS, POPSIZE, SAVEPOPS, rep, infos, seed, samp):

    nloc = loci.nloc
    pop = sim.Population(
        size        = POPSIZE, 
        loci        = loci.get_loci(),
        infoFields  = infos
    )
    ng = len(str(GENERATIONS))
    np = len(str(POPSIZE))
    ns = len(str(seed))
    sm = len(str(samp - 1))
    NG = "0"+str(ng)
    NP = "0"+str(np)
    NS = "0"+str(ns)
    SM = "0"+str(sm)
    NE = str(np + 2)

    # Init ops
    # ==========================================================================
    inits = [
        sim.Stat(effectiveSize=range(nloc), vars='Ne_temporal_base'),
        sim.InitSex(),                      # initialize sex for the populations
        sim.InitInfo(0, infoFields = infos) # set information fields to 0
    ]
    # initialize genotypes for each locus separately.
    inits += [sim.InitGenotype(prop = [0.5, 0.5])]


    # Pre Ops
    # ==========================================================================
    prelist = [
        sim.Stat(numOfMales = True, popSize = True),
        sim.TerminateIf('numOfMales == 0 or numOfMales == popSize'),
        sim.SNPMutator(u = loci.mu, v = loci.mu),
        sim.IdTagger(),
        ]

    # Post Ops
    # ==========================================================================
    # The stats operators must be raw strings. Since it would be advantageous
    # later to be able to add in different stats, it's better to join a bunch of
    # small strings than write one long one.

    # The stats must have single quotes around them.
    head, foot  = r"'", r"\n'"
    popsize     = r"Pop Size: {:"+NP+"d}"
    males       = r"Males: {:"+NP+"d}"
    Ne          = r"Ne: {n[0]:"+NE+".1f} ({n[1]:"+NE+".1f} - {n[2]:"+NE+".1f})"
    het         = r"Het: " + r"{:.2f} "*nloc
    generations = r"Gen: {:"+NG+"d}"
    reps        = r"Rep: {:2d}"

    # Joining the statistics together with pipes.
    stats = " | ".join([head, popsize, males, generations, reps, Ne, foot])

    # The string for the evaluation of the stats.
    stateval = ".format("
    stateval += "popSize, "
    stateval += "numOfMales, "
    stateval += "gen, "
    stateval += "rep, "
    stateval += "n = Ne_waples89_P1" # This returns a list of estimate and SE
    stateval += ")"

    # Stat and PyEval are both classes, so they can be put into variables. These
    # will be evaluated as the program runs. 
    statargs = sim.Stat(
        popSize = True, 
        numOfMales = True, 
        # heteroFreq = range(nloc), 
        effectiveSize = range(nloc),
        step = STEPS,
        vars = 'Ne_waples89_P1'
    )
    evalargs = sim.PyEval(stats + stateval, step = STEPS)

    postlist = [
        statargs,
        evalargs#,
        # sim.PyOperator(func = reassign_parents, step = STEPS),
    ]

    finallist = []

    if SAVEPOPS is True:
        seedf     = "seed_{:" + NS + "d}"
        sexf      = "_sex_{:1.4f}"
        sexseed   = seedf.format(seed) + sexf.format(sexrate)
        outfile   = "!'"+ sexseed + "_gen_{:"+NG+"d}_rep_{:02d}.pop'"
    # Argument Compilation
    # ==========================================================================
    # The arguments to sim.evolve can be stored as a dict and be passed as 
    # **kwargs (keyword arguments).

    EVOL = dict()
    EVOL['initOps']      = inits
    EVOL['matingScheme'] = mix_mating(sexrate)
    EVOL['preOps']       = prelist
    EVOL['postOps']      = postlist
    EVOL['finalOps']     = finallist
    EVOL['gen']          = GENERATIONS

    # Simulate and Evolve
    # ==========================================================================
    sim.IdTagger().reset(1) # IdTagger must be reset before evolving.
    simu = sim.Simulator(pop, rep = rep)

    # Saving initial simulations
    repcount = 0
    for pop in simu.populations():
        pop.SavePopulation(output = outfile + ".format(0, " + repcount + ")")
        repcount += 1

    # Print a description of the evolution process for debugging
    print(sim.describeEvolProcess(**EVOL))
    evol = simu.evolve(**EVOL)


    # Clean up populations that haven't evolved
    # ==========================================================================
    #
    # When using a mixed-mating scheme with very little sex, it's possible that
    # the population will run out of either sex before the number of generations
    # is complete. In this case, we will need to run out the rest of the
    # populations with clonal reproduction only.
    #
    # The return value for the evolution of a simulator is the number of
    # generations that were evolved for each replicate. Here, we are grabbing
    # the populations that didn't evolve to the correct number of generations.
    pops = [x for x, y in zip(simu.populations(), evol) if y < GENERATIONS]
    fin  = [x for x, y in zip(simu.populations(), evol) if y >= GENERATIONS]
    gens = [x for x in evol if x < GENERATIONS]

    foutfile = "!'"+ sexseed + "_gen_{:"+NG+"d}_rep_{:02d}'"
    saveSampleReps(simu, evol, GENERATIONS, samp, foutfile, SM)

    # 
    # Because we only want to finish out the evolution here, we put a cap on the
    # number of generations needed to evolve.
    EVOL['preOps']        = prelist[:1] + prelist[2:]
    EVOL['postOps']      += [sim.TerminateIf("gen >= "+str(GENERATIONS))]
    EVOL['matingScheme']  = sim.HomoMating(
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
    resim = sim.Simulator(pops)
    evol2 = resim.evolve(**EVOL)

    saveSampleReps(resim, evol2, GENERATIONS, samp, foutfile, SM)

# ==============================================================================
# MAIN PROGRAM
# ==============================================================================

if __name__ == '__main__':
    pars = parser.parse_args()
    
    # This variable is used to set up information fields that are collected at
    # mating.
    infos = [
        'clone_proj', # the average number of clonal events for each parent
        'sex_proj',   # average number of sexual events for each parent
        'ind_id',     # ID for each individual per generation
        'mother_id',  # ID of the mother (0 if no mother) [SIMUPOP PARAM]
        'father_id',  # ID of father (0 if no father)     [SIMUPOP PARAM]
        'tsmrsr'      # time to most recent sexual reproduction
    ]
    if not os.path.isdir(pars.outfile):
        os.mkdir(pars.outfile)
    cwd = os.getcwd()
    os.chdir(pars.outfile)

    for s in range(pars.nseed):
        loci = ra.snps(nloc = pars.nloc, nchrom = pars.nchrom, mu = pars.murate)
        for i in pars.sexrate:
            sim_partial_clone(loci, i, pars.STEPS, pars.GENERATIONS, \
                pars.POPSIZE, True, pars.rep, infos, s, pars.sample)
    args_to_file(pars)

    os.chdir(cwd)
