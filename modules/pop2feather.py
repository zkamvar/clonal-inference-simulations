#!/usr/bin/env python3.4
import numpy as np
import pandas as pd
import sys, os, re
import feather
from simuPOP import loadPopulation

def trim_lociNames(pop):
    return([l[:-2] for l in pop.lociNames()])

def get_generation(pops, gen = 10):
    rs   = r"^.+?gen_" + str(gen) + ".+?\.pop$"
    pops = [re.findall(rs, x) for x in pops]
    pops = [x[0] for x in pops if len(x) > 0]
    return(pops)

def pop2dictlist(pop, popname = None):
    '''
    Convert pop data to a dataframe

    Parameters:
        pop a simuPOP population

    Output:
        a list of dictionaries containing data for each individual in a population

    Example:
    
    import zk_utils as zu
    import simuPOP as sim
    import pandas as pd
    import os, re
    
    # Get a list of populations that are at the final generation
    pops = ["testpop/" + x for x in os.listdir("testpop")]
    pops = zu.get_generation(pops, gen = 10)

    # Load the first population and create a pandas data frame
    pop  = sim.loadPopulation(pops[0])
    dl   = zu.pop2dictlist(pop, popname = pops[0])
    cols = list(pop.infoFields()) + ["sex", "pop"] + zu.trim_lociNames(pop)
    df   = pd.DataFrame(dl, columns = cols)

    '''
    nloc    = pop.totNumLoci()
    # The locus names need to be stripped of the repeat lengths
    lnames  = trim_lociNames(pop)
    infos   = pop.infoFields()
    ploid   = pop.ploidy()
    outlist = []
    for ind in pop.individuals():
        indict = {}
        indict["sex"] = ind.sex()
        if popname is not None:
            indict["pop"] = popname
        # Information fields for each individual
        for i in infos:
            indict[i] = ind.info(i)
        # locus data
        for x in range(nloc):
            # list comprehension of alleles for each homologous chromosome
            geno = "/".join([str(ind.allele(x, p) + 1) for p in range(ploid)])
            indict[lnames[x]] = geno
        outlist.append(indict)
    return(outlist)

def pops2df(pops):
    dl = []
    for p in pops:
        pop = loadPopulation(p)
        dl += pop2dictlist(pop, p)
    cols = list(pop.infoFields()) + ["sex", "pop"] + trim_lociNames(pop)
    return(pd.DataFrame(dl, columns = cols))

if __name__ == '__main__':
	if len(sys.argv) < 2:
		sys.exit()
	d = sys.argv[1]
	if not os.path.isdir(d):
		sys.exit()
	
	pops = [d + "/" + x for x in os.listdir(d)]
	if len(sys.argv) > 2:
		pops = get_generation(pops, sys.argv[2])
	df = pops2df(pops)
	feather.write_dataframe(df, d+"/"+d+".feather")
