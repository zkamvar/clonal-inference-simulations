#!/usr/bin/env python3.4

import os
import simuOpt
import multiprocessing
cpu = int(multiprocessing.cpu_count())
simuOpt.setOptions(optimized = True, 
	gui = False, 
	debug = 'DBG_WARNING',
	alleleType = 'long', 
	# quiet = True, 
	numThreads = 4)
import simuPOP as sim
import random_alleles as ra
from simuPOP.utils import export
from simuPOP.utils import saveCSV
from simuPOP.sampling import drawRandomSample
from simuPOP.sampling import drawRandomSamples 

options = [
    {'name':'POPSIZE',
     'default':1000,
     'label':'Initial Population Size',
     'type': 'integer',
     'description':'This will set the effective population size to n individuals.',
     'validator': 'POPSIZE > 0',
     },
    {'name':'nloc',
     'default': 10,
     'label':'Number of loci',
     'description':'Defines the number of loci. It is important to note that putting in the number n will put n loci on one chromosome while puting 1 n times will put n loci on each of n chromosomes.',
     'type':'integer',
     'validator':simuOpt.valueGT(0),
    },
    {'name':'outfile',
     'default':'foo',
     'label':'Please name your output files (Letters and Numbers only)',
     'description':'gives a common name across your output files. Note that this program outputs 1+32n the number of replicates',
     'type':'string',
     'validator':'True',
    },
    {'name':'GENERATIONS',
     'default':10001,
     'label':'Number of Generations',
     'type': 'integer',
     'validator': 'GENERATIONS > 0',
    },
    {'name':'STEPS',
     'default':1000,
     'label':'Steps at which to evaluate progress of mating',
     'type': 'integer',
     'validator': 'STEPS > 0',
    },
    # {'name':'SAVEPOPS',
    #  'default': False,
    #  'label':'Should populations be saved at the number of steps?',
    #  'type': 'boolean',
    #  'validator':simuOpt.valueTrueFalse(),
    # },
    {'name':'sexrate',
     'default':0,
     'label':'Percentage of sexual reproduction.',
     'description':'This defines the coefficient of sexual reproduction that the population will undergo per generation. Valid input is any integer between 0 and 100.',
     'type':'number',
     'validator':'100 >= sexrate >= 0',
    },
    {'name':'murate',
     'default':[1e-5]*10,
     'label':'mutation rate',
     'description':'This defines the coefficient of sexual reproduction that the population will undergo per generation. Valid input is any integer between 0 and 100.',
     'type':'numbers',
     'validator': simuOpt.valueListOf(simuOpt.valueBetween(0, 1)),
    },

    {'name':'amin',
     'default':6,
     'label':'min number of alleles',
     'type':'integer',
     'validator':'amin > 0',
    },
    {'name':'amax',
     'default':10,
     'label':'max number of alleles',
     'type':'integer',
     'validator':'amax > amin',
    },
    {'name':'rep',
     'default':2,
     'label':'Number of replicates populations',
     'type': 'integer',
     'validator': 'rep > 0',
    }
]


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

def sim_partial_clone(loci, sexrate, STEPS, GENERATIONS, POPSIZE, SAVEPOPS, rep):
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
	nloc = loci.nloc()
	pop = sim.Population(
	    size        = POPSIZE, 
	    loci        = [1]*nloc, 
	    lociNames   = loci.get_locus_names(), 
	    alleleNames = loci.get_allele_names(),
	    infoFields  = infos
	)

	# Init ops
	# ==========================================================================
	freqs = loci.get_frequencies()
	inits = [
		sim.InitSex(),                      # initialize sex for the populations
		sim.InitInfo(0, infoFields = infos) # set information fields to 0
	]
	# initialize genotypes for each locus separately.
	inits += [sim.InitGenotype(freq = freqs[i], loci = i) for i in range(nloc)]


	# Mating Schemes
	# ==========================================================================
	# There are options here are things to do during mating:
	# 1. Tag the parents (which fills mother_idx and father_idx)
	# 2. Count the reproductive events

	rand_mate = sim.RandomMating(
	    numOffspring=(sim.UNIFORM_DISTRIBUTION, 1, 3),
	    ops = [
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
	mate_scheme = sim.HeteroMating([rand_mate, clone_mate])


	# Post Ops
	# ==========================================================================
	# The stats operators must be raw strings. Since it would be advantageous
	# later to be able to add in different stats, it's better to join a bunch of
	# small strings than write one long one.

	# The stats must have single quotes around them.
	head, foot  = r"'", r"\n'"
	popsize     = r"Pop Size: %d"
	males       = r"Males: %d"
	het         = r"Het: " + r"%.2f "*nloc
	generations = r"Gen: %d"

	# Joining the statistics together with pipes.
	stats = " | ".join([head, popsize, males, het, generations, foot])

	# Heterozygosity must be evaluate for each locus. This is a quick and dirty
	# method of acheiving display of heterozygosity at each locus. 
	locrange = map(str, range(nloc))
	lochet = '], heteroFreq['.join(locrange)

	# The string for the evaluation of the stats.
	stateval = " % (popSize, numOfMales, heteroFreq["+lochet+"], gen)"

	# Stat and PyEval are both classes, so they can be put into variables. These
	# will be evaluated as the program runs. 
	statargs = sim.Stat(
		popSize = True, 
		numOfMales = True, 
		heteroFreq = range(nloc), 
		step = STEPS
	)
	evalargs = sim.PyEval(stats + stateval, step = STEPS)

	postlist = [
		statargs,
		evalargs#,
		# sim.PyOperator(func = reassign_parents, step = STEPS),
	]

	if SAVEPOPS is True:
	    outfile = "!'gen_%d.pop' % (gen)"
	    postlist += [sim.SavePopulation(output = outfile, step = STEPS)]


	# Simulate and Evolve
	# ==========================================================================	
	sim.IdTagger().reset(1) # IdTagger must be reset before evolving.
	simu = sim.Simulator(pop, rep = rep)
	simu.evolve(
	    initOps = inits,
	    matingScheme = mate_scheme,
	    preOps = [
	        sim.StepwiseMutator(rates = loci.get_mu(), loci = range(nloc)),
	        sim.IdTagger(),
	        ],
	    postOps = postlist,
	    gen = GENERATIONS
	)


if __name__ == '__main__':
    pars = simuOpt.Params(options, 'OPTIONS')
    if not pars.getParam():
        sys.exit(0)

    # First, generate the 
    loci = generate_loci(pars.nloc, pars.murate, pars.amax, pars.amin)

    sim_partial_clone(loci, pars.sexrate, pars.STEPS, pars.GENERATIONS, \
    	pars.POPSIZE, False, pars.rep)

    pars.saveConfig("config.cfg")
