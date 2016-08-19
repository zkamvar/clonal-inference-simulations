#!/usr/bin/env python3.4

# Useful for importing any python modules I have in the "modules" directory.
# http://stackoverflow.com/a/35259170/2752888
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'modules'))
import random_alleles as ra
import argparse
import inspect

def convert_arg_line_to_args(arg_line):
    for arg in arg_line.split():
        if not arg.strip():
            continue
        if arg[0] == '#':
            break
        yield arg


def args_to_file(args):
    f = open(args.cfg, "w")
    for arg, value in vars(args).items():
        if not isinstance(value, (list)):
            value = [value]
        val = '{} '*len(value)
        f.write("--{} ".format(arg) + val.format(*value) + "\n")
    f.close


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
    help = "Set the number of unlinked loci.",
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
    nargs = "+",
    help = "mutation rate per locus",
    default = [1e-05]*10
    )
parser.add_argument(
    "--amin", 
    type = int,
    help = "Minimum number of alleles per locus",
    default = 6
    )
parser.add_argument(
    "--amax", 
    type = int,
    help = "Maximum number of alleles per locus",
    default = 10
    )
parser.add_argument(
    "--rep", 
    type = int,
    help = "Number of replicates per population",
    default = 10
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


def generate_loci(nloc, murate, amax, amin):
    # Initializing the allele frequencies for each locus.
    #
    # Each element in 'alleles' is a class of zk_locus, which is then fed into
    # zk_loci, a class that contains zk_locus and provides access to names, 
    # frequencies, and mutation rates.
    alleles = []
    for x in range(nloc):
        alleles += [ra.zk_locus(mu = murate[x], amax = amax, amin = amin)]
    loci = ra.zk_loci(alleles)
    return(loci)

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


def sim_partial_clone(loci, sexrate, STEPS, GENERATIONS, POPSIZE, SAVEPOPS, rep, infos, seed):

    nloc = loci.nloc()
    pop = sim.Population(
        size        = POPSIZE, 
        loci        = [1]*nloc, 
        lociNames   = loci.get_locus_names(), 
        alleleNames = loci.get_allele_names(),
        infoFields  = infos
    )
    ng = len(str(GENERATIONS))
    np = len(str(POPSIZE))
    NG = "0"+str(ng)
    NP = "0"+str(np)
    NE = str(np + 2)

    # Init ops
    # ==========================================================================
    freqs = loci.get_frequencies()
    inits = [
        sim.Stat(effectiveSize=range(nloc), vars='Ne_temporal_base'),
        sim.InitSex(),                      # initialize sex for the populations
        sim.InitInfo(0, infoFields = infos) # set information fields to 0
    ]
    # initialize genotypes for each locus separately.
    inits += [sim.InitGenotype(freq = freqs[i], loci = i) for i in range(nloc)]


    # Pre Ops
    # ==========================================================================
    prelist = [
        sim.Stat(numOfMales = True, popSize = True),
        sim.TerminateIf('numOfMales == 0 or numOfMales == popSize'),
        sim.StepwiseMutator(rates = loci.get_mu(), loci = range(nloc)),
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
        sexf      = "seed_{:2d}_sex_{:1.4f}".format(seed, sexrate)
        outfile   = "!'"+sexf+"_gen_{:"+NG+"d}_rep_{:02d}.pop'.format(gen, rep)"
        postlist += [sim.SavePopulation(output = outfile, step = STEPS)]
        # finallist += [sim.SavePopulation(output = outfile)]

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
    gens = [x for x in evol if x < GENERATIONS]
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
    resim.evolve(**EVOL)

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
        loci = generate_loci(pars.nloc, pars.murate, pars.amax, pars.amin)
        for i in pars.sexrate:
            sim_partial_clone(loci, i, pars.STEPS, pars.GENERATIONS, \
                pars.POPSIZE, True, pars.rep, infos, s)
    args_to_file(pars)

    os.chdir(cwd)
