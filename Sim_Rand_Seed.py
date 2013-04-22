#!/usr/bin/env python

# This script was written by Zhian Kamvar

import simuOpt, os, sys, datetime, commands, re, exceptions, multiprocessing, random, pickle
from simuOpt import setOptions
cpu = int(multiprocessing.cpu_count())
setOptions(optimized=True, gui=True, debug='DBG_WARNING',alleleType='long', quiet=True, numThreads=cpu)
import simuPOP as sim
from simuPOP import utils
from simuPOP.sampling import drawRandomSample
#------------------------------------------------------------------------------#
# setting the options for the script:
# S0 = Pop size
# Sam0 = Sample Size
# G0 = Number of Generations
# C0 = Percent of clonal reproduction
# R0 = Number of Replicates
# L0 = Number of loci
# N0 = Name of the fileset
#------------------------------------------------------------------------------#
options = [
    {'name':'S0',
     'default':10000,
     'label':'Initial Population Size',
     'type': 'integer',
     'description':'This will set the effective population size to n individuals.',
     'validator': 'S0 > 0',
     },
    {'name':'Sam0',
     'default': (10, 25, 50, 100),
     'label':'Sample Population Size',
     'type': 'integers',
     'description':'Chooses random samples from the population.',
     'validator': simuOpt.valueListOf(simuOpt.valueGT(0)),
    },
    {'name':'G0',
     'default':10000,
     'label':'Number of Generations',
     'type': 'integer',
     'validator': 'G0 > 0',
    },
    {'name':'C0',
     'default':100,
     'label':'Percentge of clonal reproduction.',
     'description':'This defines the coefficient of clonal reproduction that the population will undergo per generation. Valid input is any integer between 0 and 100.',
     'type':'number',
     'validator':'100 >= C0 >= 0',
    },
    {'name':'R0',
     'default':1,
     'label':'Number of Replicate SAMPLES you wish to take. If you have 4 numbers in Sam0 and you want 10 replicates, you will get 40 *.dat files',
     'type': 'integer',
     'validator': 'R0 > 0',
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
#------------------------------------------------------------------------------#
# Naming strategy:
# Given that the number of replicates and generations will change, our old
# strategy of naming them foo_0_0.csv was uninformative unless you knew what
# the program did.
# 
# New Naming Strategy:
# foo_rep_00_gen_1.00k.dat
# foo_rep_00_gen_1000.pop
# This is more informative and should convey what is needed. The function below 
# will determine the suffix needed for the .csv files depending on the number
# of generations.
#------------------------------------------------------------------------------#
def saveFStat(pop, output='', maxAllele=0, loci=[], shift=1,
    combine=None):
    '''
    '''
    if output != '':
        file = output
    else:
        raise exceptions.ValueError, "Please specify output"
    # open file
    try:
        f = open(file, "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + file + " to write."
    #
    # file is opened.
    np = pop.numSubPop()
    if np > 200:
        print "Warning: Current version (2.93) of FSTAT can not handle more than 200 samples"
    if loci == []:
        loci = range(pop.totNumLoci())
    nl = len(loci)
    if nl > 100:
        print "Warning: Current version (2.93) of FSTAT can not handle more than 100 loci"
    if maxAllele != 0:
        nu = maxAllele
    else:
        nu = max(pop.genotype()) + shift
    if nu > 999:
        print "Warning: Current version (2.93) of FSTAT can not handle more than 999 alleles at each locus"
        print "If you used simuPOP_la library, you can specify maxAllele in population constructure"
    if nu < 10:
        nd = 1
    elif nu < 100:
        nd = 2
    elif nu < 1000:
        nd = 3
    else: # FSTAT can not handle this now. how many digits?
        nd = len(str(nu))
    # write the first line
    f.write( '%d %d %d %d\n' % (np, nl, nu, nd) )
    # following lines with loci name.
    for loc in loci:
        #
        f.write( pop.locusName(loc) +"\n");
    gs = pop.totNumLoci()
    for sp in range(0, pop.numSubPop()):
        # genotype of subpopulation sp, individuals are
        # rearranged in perfect order
        gt = pop.genotype(sp)
        for ind in range(0, pop.subPopSize(sp)):
            f.write("   %d  " % (sp+1))
            p1 = 2*gs*ind          # begining of first hemo copy
            p2 = 2*gs*ind + gs     # second
            for al in loci: # allele
                # change from 0 based allele to 1 based allele
                if combine is None:
                    ale1 = gt[p1+al] + shift
                    ale2 = gt[p2+al] + shift
                    f.write('%%0%dd%%0%dd ' % (nd, nd) % (ale1, ale2))
                else:
                    f.write('%%0%dd' % nd % combine([gt[p1+al], gt[p2+al]]))
            f.write( "\n")
    f.close()

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
#------------------------------------------------------------------------------#
# This is another object that has to do with naming. It will place a zero in 
# front of any number that is less than the final number of generations for 
# ease of listing.
#------------------------------------------------------------------------------#
def regul(m,G0):
    if (m < G0):
        m = "0%.1f" % (float(m)/upp)
        return m
    elif (m == G0):
        m = "%.1f" % (float(m)/upp)
        return m
    else:
        return m
#------------------------------------------------------------------------------#
# Returns all possible loci names
#------------------------------------------------------------------------------#
def loci_list(all_loci):
    for i in all_loci:
        all_loci[i] = "Locus_"+str(i+1)
    return all_loci
#------------------------------------------------------------------------------#
# clone simulator in case the population goes strictly clonal.
# newG0 is the number of generations left for the population to go.
# This will also create a population marking the point at which all individuals
# are of the same sex (or mating type). 
#------------------------------------------------------------------------------#
def goToClone(pops, L0, newG0, G0, N0, R0):   
    simu = sim.Simulator(pops)
    print("GOING TO CLONAL REPRODUCTION")
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

def newRandomizer():
    num_of_alleles = random.randint(5,10)
    new = list(range(10))
    for i in range(10):
        if i < num_of_alleles:
            new[i] = random.randint(1,10)
        else:
            new[i] = 0
    listsum = float(sum(new))
    for i in range(num_of_alleles):
        new[i] = new[i]/listsum
    print(new)
    return(new)

def randLoci(loci):
    loclist = (list(range(loci)))
    for i in range(loci):
       loclist[i] = newRandomizer()
    return(loclist) 
#------------------------------------------------------------------------------#
# simulation object
#------------------------------------------------------------------------------#


def partialCloneFromBurnin(Sam0, S0, G0, C0, R0, N0, L0):
    repnum = N0.split("_")[5]
    pop = sim.loadPopulation("./BURNIN_pop_"+str(S0)+"_rep_"+repnum+".pop")
    sample = drawRandomSample(pop, sizes=Sam0)
    commands.getoutput('mkdir '+str(N0))
    print("Percent Random Mating: "+str(100-C0)+"\nPercent Clone Mating: "+str(C0))
		
#------------------------------------------------------------------------------#
# Everything above initiates the population, sample size, and creates the
# directory for all of the data.
#------------------------------------------------------------------------------#
    simu = sim.Simulator(pop)
    gens = simu.evolve(
        # This is supposed to give proportional weights to each scheme, so if 
        # C0 = 25, then rand = 75 and the number of individuals (N) that will 
        # be affected by this will be (25/100)*N and (75/100)*N, respectively.  
        matingScheme = sim.HeteroMating([
            # The following two mating schemes perform clonal mating and they
            # choose an equal number of male and female individuals for mating,
            # so that the simulation will never run out. 
            sim.RandomSelection(subPops = 0, weight = C0),
            # RandomMating is set to assign sex randomly in equal proportions by
            # default. 
            #sim.IfElse('numOfMales!=0 and numOfFemales !=0',
            #    ifOps = sim.RandomMating(subPops = 0, weight = (100-C0)),
            #    elseOps = sim.RandomSelection(subPops = 0, weight = (100-C0))
            #),
            sim.RandomMating(subPops=0, weight=(100-C0))
        ]),
        preOps=[
            # Since this simulation will still have the random mating as part
            # of the mating scheme, it will return an error if it runs out of 
            # either sex. This terminator will end the simulation if this
            # happens and functions below will start the population as clones
            sim.Stat(numOfMales=True),
            sim.TerminateIf(condition='numOfMales == 0 or numOfFemales == 0'),
            sim.StepwiseMutator(rates=1e-5, loci=L0)
        ],
        postOps=[
            # This outputs more informative statistics to the screen. 
            sim.Stat(popSize=True, numOfMales=True, heteroFreq=L0, step=G0/10),
            sim.PyEval(r"'Pop Size: %d | Number of Males: %d | Heterozygosity: \
%.2f | Gen: %d\n' % (popSize, numOfMales, heteroFreq[1], gen-("+str(G0)+"/10))", step=G0/10),
            
            #sim.Stat(LD=[0, 1], step=G0/10),
            #sim.PyEval(r"'%.2f\n' % LD[0][1]", step=G0/10),
#------------------------------------------------------------------------------#
# In the postOps, the two functions above are useless at the moment.
# Below, the sim.SavePopulation saves each population in increments of 1/10 the
# number of generations
#------------------------------------------------------------------------------#
            sim.SavePopulation(output="!'"+str(N0)+"/"+str(N0)+
            "_gen_%d.pop' % (gen-("+str(G0)+"/10))", step=G0/10),
        ],
        finalOps=[
            sim.SavePopulation(output="!'"+str(N0)+"/"+str(N0)+
            "_gen_%d.pop' % (gen-("+str(G0)+"/10))", step=G0/10),
        ],
        gen = G0
    )
    # Taking the return value (gens) from evolve and putting the populations
    # that did not complete the full allotted generations here.
    pops = [x for x, y in zip(simu.populations(), gens)]#if y < G0]
    #print(pops)
    #print(simu.numRep())
    #print(zip(simu.populations()))
    # Creating a list of unfinished generation times.
    unfinished = [x for x, y in zip(gens, gens)]#if y < G0]
    print("Unfinished: ")
    print(unfinished)
    # making a loop over the number of unfinished populations and iterating a 
    # clonal mating scheme over them until they finish.
    for i in range(len(unfinished)):
        if G0 == unfinished[i]:
            continue
        print("Rep "+str(i+1)+" needs more work!\n")
        newG0 = G0 - unfinished[i]
        goToClone(pops[i], L0, newG0, G0, N0, i+1)
#------------------------------------------------------------------------------#
# After all the simulations are done, a sample of the final population is drawn
# and then all of the previous saved populations are sampled from using the
# sim.loadPopulation function and all are saved as csv files
#------------------------------------------------------------------------------#
    for m in range(0,G0+1,(G0)/10):
        ppop = sim.loadPopulation(str(N0)+"/"+str(N0)+"_gen_%d.pop" % m)
        for k in range(R0):
            for a in range(0,len(Sam0)):
                sam = drawRandomSample(ppop, sizes=Sam0[a])
                b = regul(m,G0)
#                print(Sam0[a])
                outfile = str(N0)+"/"+str(N0)+"_sr_%02d" % (k+1)+"_gen_"+str(b)+str(nam)+"_sam_"+str(Sam0[a])+".dat"
                saveFStat(sam, outfile)
    return



#------------------------------------------------------------------------------#
# Finally, we get to the actual script itself. The first action is to make sure 
# that the program is not called from within another script (or something like
# that) with the if __name__ = '__main__' thing. It then checks if the filename
# you want exists.
#------------------------------------------------------------------------------#
if __name__ == '__main__':
    pars = simuOpt.Params(options, 'OPTIONS')
    if not pars.getParam():
        sys.exit(0)
    if os.path.exists(r'./'+pars.N0):
        print (
            "\n\n######################\n\n"
            "ERROR:\nFILENAME "+pars.N0+" EXISTS.\nPLEASE CHOOSE ANOTHER.\n\n"
            "######################\n\n"
        )
        sys.exit(0)
#------------------------------------------------------------------------------#
# This is where the configuration file and "LOG" files are written/initiated
#------------------------------------------------------------------------------#
    pars.saveConfig(pars.N0+".cfg")
    (nam, upp) = testnum(pars.G0)
    tim = str(datetime.datetime.now())
    with open('./LOG_'+str(pars.N0)+'.txt', 'w') as f:
        f.write(
            "\n#------------------------------------------------------------------------------#\n"
            "#\tRun started on "+tim+"\n"
            )
        f.close()
    #loclist = randLoci(len(pars.L0))
    repnum = pars.N0.split("_")[5]
    loclist = pickle.load(open("BURNIN_pop_"+str(pars.S0)+"_rep_"+str(repnum)+".p", "rb"))
#------------------------------------------------------------------------------#
# The running of the simulation with all the parameters is given by this line.
#------------------------------------------------------------------------------#
    partialCloneFromBurnin(pars.Sam0, pars.S0, pars.G0, pars.C0, pars.R0, pars.N0, pars.L0)
#------------------------------------------------------------------------------#
# It finishes by writing out the rest of the "LOG" file with all the parameters.
#------------------------------------------------------------------------------#
    tim = str(datetime.datetime.now())
    with open('./LOG_'+str(pars.N0)+'.txt', 'a') as f:
        count = commands.getoutput('ls '+str(pars.N0)+' | grep -c '+str(pars.N0))
        f.write(
            "#\tRun finished on "+tim+"\n#\n"
            "#\tThe following simuPOP parameters were used in this simulation:\n"
            "#\tPopulation size: "+str(pars.S0)+"\n"
            "#\tSample size: "+str(pars.Sam0)+"\n"
            "#\tNumber of Generations: "+str(pars.G0)+"\n"
            "#\tPercent clonal Reproduction: "+str(pars.C0)+"\n"
            "#\tNumber of Sampling Replicates: "+str(pars.R0)+"\n"
            "#\tNumber of Loci: "+str(pars.L0)+"\n"
            "#\tFileset name: "+str(pars.N0)+"\n"
            "#\tTotal number of files produced: "+str(count)+"\n"
            "#\tAllele Frequencies Per locus: \n"
            )
        for i in loclist:
            counts = 0
            f.write("#\t")
            for j in i:
                counts += 1
                if counts < 6:
                   f.write(str(j)+" ")
                else:
                    counts = 0
                    f.write("\n#\t"+str(j)+" ")
                            
            f.write("\n\n")
        f.write(
            "#\n"
            "#\tProcessors Used: "+str(cpu)+"\n"
            "#------------------------------------------------------------------------------#\n"
            )
        f.close()
