#!/usr/bin/env python

# This script was written by Zhian Kamvar

import simuOpt, os, sys, datetime, commands, re, exceptions, multiprocessing, random, pickle
from simuOpt import setOptions
cpu = int(multiprocessing.cpu_count())
setOptions(optimized=True, gui=True, debug='DBG_WARNING',alleleType='long', quiet=True, numThreads=cpu)
import simuPOP as sim
from simuPOP import utils
from simuPOP.sampling import drawRandomSample
from random_alleles import *
options = [
    {'name':'S0',
     'default':10000,
     'label':'Initial Population Size',
     'type': 'integer',
     'description':'This will set the effective population size to n individuals.',
     'validator': 'S0 > 0',
     },
    {'name':'B0',
     'default':1000,
     'label':'Number of Generations to Burn in',
     'type': 'integer',
     'validator': 'B0 > 0',
    },
    {'name':'L0',
     'default':([1]*10),
     'label':'Loci information',
     'description':'Defines the number of loci. It is important to note that putting in the number n will put n loci on one chromosome while puting 1 n times will put n loci on each of n chromosomes.',
     'type':'integers',
     'validator':simuOpt.valueListOf(simuOpt.valueGT(0)),
    },
    {'name':'N0',
     'default':'foo',
     'label':'Please name your output files (Letters and Numbers only)',
     'description':'gives a common name accross your output files. Note that this program outputs 1+32n the number of replicates',
     'type':'string',
     'validator':'True',
    },
]

# shortening number strings for file names. 

def testnum(num):
    if len(str(num)) > 3:
        nam = "k"
        upp = 1000
        if len(str(num)) > 6:
            nam = "mil"
            upp = 1000000
            if len(str(num)) > 9:
                nam = "bil"
                upp = 1000000000
                return(nam, upp)
            else:
                return(nam, upp)
        else:
            return(nam, upp)
    else:
        nam = ""
        upp = 1
        return(nam, upp)

# 

def regul(m,G0):
    if (m < G0):
        m = "0%.1f" % (float(m)/upp)
        return m
    elif (m == G0):
        m = "%.1f" % (float(m)/upp)
        return m
    else:
        return m

def loci_list(all_loci):
    for i in all_loci:
        all_loci[i] = "Locus_"+str(i+1)
    return all_loci

def newRandomizer():
    num_of_alleles = random.randint(6,10)
    new = list(range(10))
    for i in range(10):
        if i < num_of_alleles:
            new[i] = random.randint(1,10)
        else:
            new[i] = 0
    listsum = float(sum(new))
    for i in range(num_of_alleles):
        new[i] = new[i]/listsum
    #print(new)
    return(new)

def randLoci(loci):
    loclist = (list(range(loci)))
    for i in range(loci):
       loclist[i] = newRandomizer()
    return(loclist) 

def burnin(S0, B0, L0, loclist, N0):
    all_loci = list(range(len(L0)))
    allNames=('1','2','3','4','5','6','7','8','9','10')
    pop = sim.Population(
        size=S0,
        loci = L0,#[int(s) for s in L0.split()],break loop python
        lociNames = (loci_list(all_loci)),
        alleleNames=('1','2','3','4','5','6','7','8','9','10'),
        ploidy=2,)
#    sample = drawRandomSample(pop, sizes=Sam0)
#    commands.getoutput('mkdir '+str(N0))
    #alProps = randomDistPerRep(allNames, N0, R0)
    loclist = loclist
    initlist = list()
    initlist.append(sim.InitSex(maleFreq = 0.5))
    for i in range(len(loclist)):
        initlist.append(sim.InitGenotype(freq = loclist[i], loci = i))
    simu = sim.Simulator(pop)
    simu.evolve(
        initOps= initlist,
        matingScheme = sim.RandomMating(),
        finalOps = sim.SavePopulation(output="!'"+str(N0)+".pop'"),
        gen = B0
    )

if __name__ == '__main__':
    pars = simuOpt.Params(options, 'OPTIONS')
    if not pars.getParam():
        sys.exit(0)
    pars.saveConfig(str(pars.N0)+".cfg")
    (nam, upp) = testnum(pars.B0)
    tim = str(datetime.datetime.now())
    with open('./LOG_'+str(pars.N0)+'.txt', 'w') as f:
        f.write(
            "\n#------------------------------------------------------------------------------#\n"
            "#\tBurnin started on "+tim+"\n"
            )
        f.close()
    loclist = get_allele_probabilities(len(pars.L0), 10)
    pickle.dump(loclist, open(str(pars.N0)+".p", "wb"))
#------------------------------------------------------------------------------#
# The running of the simulation with all the parameters is given by this line.
#------------------------------------------------------------------------------#
    burnin(pars.S0, pars.B0, pars.L0, loclist, pars.N0)
#------------------------------------------------------------------------------#
# It finishes by writing out the rest of the "LOG" file with all the parameters.
#------------------------------------------------------------------------------#
    tim = str(datetime.datetime.now())
    with open('./LOG_'+str(pars.N0)+'.txt', 'a') as f:
        f.write(
            "#\tRun finished on "+tim+"\n#\n"
            "#\tThe following simuPOP parameters were used in this BURNIN:\n"
            "#\tPopulation size: "+str(pars.S0)+"\n"
            "#\tNumber of Burn in Generations: "+str(pars.B0)+"\n"
            "#\tNumber of Loci: "+str(pars.L0)+"\n"
            "#\tAllele Frequencies Per locus: \n"
            )
        for i in loclist:
            count = 0
            f.write("#\t")
            for j in i:
                count += 1
                if count < 6:
                   f.write(str(j)+" ")
                else:
                    count = 0
                    f.write("\n#\t"+str(j)+" ")
                            
            f.write("\n\n")
        f.write(
            "#\n"
            "#\tProcessors Used: "+str(cpu)+"\n"
            "#------------------------------------------------------------------------------#\n"
            )
        f.close()
