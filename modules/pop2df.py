#!/usr/bin/env python3

import sys, os, re
# Setting up options
import simuOpt
simuOpt.setOptions(optimized = False, 
    gui = False, 
    debug = 'DBG_WARNING',
    alleleType = 'long', 
    quiet = False, 
    numThreads = 0)
import simuPOP as sim
import numpy as np
import pandas as pd


def trim_lociNames(pop):
    '''
    Trim the last two elements of the locus names
    '''
    return([l[:-2] for l in pop.lociNames()])

def get_field(pops, regex = "gen_10"):
    '''
    Return all strings that have a specific regex

    pops = ["testpop/" + x for x in os.listdir("testpop")]

    # Get everything that has completed generation times
    pops = zu.get_field(pops, "gen_10")
    '''
    rs   = re.compile("^.+?" + str(regex) + ".*?\.pop$")
    pops = [re.findall(rs, x) for x in pops]
    pops = [x[0] for x in pops if len(x) > 0]
    return(pops)

def pop2dictlist(pop, popname = None):
    '''
    Convert pop data to a list of dictionaries

    Parameters:
        pop a simuPOP population

    Output:
        a list of dictionaries containing data for each individual in a population

    Example:
    
    import pop2df
    import simuPOP as sim
    import pandas as pd
    import os, re
    
    # Get a list of populations that are at the final generation
    pops = ["testpop/" + x for x in os.listdir("testpop")]
    pops = pop2df.get_field(pops, "gen_10")

    # Load the first population and create a pandas data frame
    pop  = sim.loadPopulation(pops[0])
    dl   = pop2df.pop2dictlist(pop, popname = pops[0])

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
    '''
    Convert a list of populations to a pandas DataFrame

    Parameters:
        pops a list of simuPOP populations

    Output:
        a pandas DataFrame

    Example:

    import pop2df
    import simuPOP as sim
    import pandas as pd
    import os, re

    # Get a list of populations that are at the final generation
    pops = ["testpop/" + x for x in os.listdir("testpop")]
    pops = pop2df.get_field(pops, "gen_10")

    df   = pop2df.pops2df(pops)
    '''
    dl     = []
    INFO   = False
    infos  = []
    lnames = []
    # print(os.getcwd())
    for p in pops:
        pop = sim.loadPopulation(p)
        dl += pop2dictlist(pop, p)
        if not INFO:
            lnames += trim_lociNames(pop)
            infos  += pop.infoFields()
            INFO   = True
    cols = infos + ["sex", "pop"] + lnames
    return(pd.DataFrame(dl, columns = cols))
