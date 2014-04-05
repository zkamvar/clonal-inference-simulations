#!/usr/bin/env python

import os, sys, getopt


def usage():
    print("\n%s Copyright (C) Zhian N. Kamvar 2014" % sys.argv[0])
    print("Software comes with no warranty.")
    print("This will generate files for simuPOP simulations.")
    print("Usage:\n\tpython %s -d '.' -r 100 -l 10 -p '10000' -s '10 25 50 100'" % sys.argv[0])
    print("Options:")
    print("\t-h:\tShow help and quit.")
    print("\t-v:\tVerbose.")
    print("\t-d: Directory to dump files.")
    print("\t-r: Number of replicate files.")
    print("\t-p: Population size(s) (in quotes, space separated)")
    print("\t-s: Sample size(s) (in quotes, space separated)")
    print("")
    sys.exit(2)



def write(N0, C0, S0, R0, Sam0, L0):
    with open("./clone_"+str(N0)+"_pop_"+str(S0)+"_rep_%02d" % R0+".cfg", 'w') as g:
        g.write(
            "[optimized]\n\n"
            "[S0] # Population Size\n"
            "S0="+str(S0)+"\n\n"
            "[Sam0] # Sample Size\n"
            "Sam0="+str(Sam0)+"\n\n"
            "[G0] # Number of Generations\n" 
            "G0=10000\n\n"
            "[C0] # Percent Clonal Reproduction\n"
            "C0="+'%.2f' % C0+"\n\n"
            "[R0] # Number of Replicates\n"
            "R0=10\n\n" 
            "[L0] # Number of Loci\n"
            "L0=[1]*"+str(L0)+"\n\n"
            "[N0] # Fileset Name\n"
            "N0=clone_"+str(N0)+"_pop_"+str(S0)+"_rep_%02d" % R0+"\n"
            )
        g.close()


def burn(N0, S0, R0, L0):
    with open("./BURNIN_pop_"+str(S0)+"_rep_%02d" % R0+".cfg", 'w') as g:
        g.write(
            "[optimized]\n\n"
            "[S0] # Population Size\n"
            "S0="+str(S0)+"\n\n"
            "[B0] # Number of Burnin Generations\n" 
            "B0=1000\n\n"
            "[L0] # Number of Loci\n"
            "L0=[1]*"+str(L0)+"\n\n"
            "[N0] # Fileset Name\n"
            "N0=BURNIN_pop_"+str(S0)+"_rep_%02d" % R0+"\n"
            )
        g.close()

def name(num):
    if num < 100:
        if num < 10:
            num = "00%.2f" % num
            return num
        else:
            num = "0%.2f" % num
            return num
    else:
        num = "%.2f" % num
        return num

        
if __name__ == '__main__':

    clones = [100, 99.99, 99.95, 99.9, 99, 95, 90, 80, 50, 0]
    replicates = 100
    loci = 10
    popsize = '10000'
    verbose = True
    sampsize = "10 25 50 100"
    directory = "."

    myopts, args = getopt.getopt(sys.argv[1:], "hvd:r:l:p:s:", ["directory=","replicates=", "loci=", "popsize=", "sampsize="])

    for opt, arg in myopts:
        if opt == "-h":
            usage()
        elif opt == "-v":
            verbose == True
        elif opt in ("-d", "--directory"):
            directory = arg
        elif opt in ("-r", "--replicates"):
            replicates == int(arg)
        elif opt in ("-l", "--loci"):
            loci == int(arg)
        elif opt in ("-p", "--popsize"):
            popsize == arg
        elif opt in ("-s", "--sampsize"):
            sampsize == arg

    sampsize = map(int, sampsize.split())
    popsize = map(int, popsize.split())

    if not os.path.exists(directory):
        os.makedirs(directory)
    here = os.getcwd()
    os.chdir(directory)

    if verbose:
        clen = len(clones)
        samlen = len(sampsize)
        poplen = len(popsize)
        cfgs = clen * replicates * poplen
        burns = poplen * replicates
        print("\nProducing "+str(cfgs + burns)+" *.cfg files in " + os.getcwd())
        print("Population Size(s) :" + str(popsize))
        print("Number of Loci: " +str(loci))
        print("Number of replicate runs: " +str(replicates))
        print("Size of samples: " +str(sampsize))
        allfiles = cfgs * samlen + sum(popsize)
        print("Number of output files that will be produced from this run: " +str(allfiles)+ "\n")

    for S0 in popsize:
        for C0 in clones:
            N0 = name(C0)
            for R0 in range(replicates):
                write(N0, C0, S0, R0, sampsize, loci)
                if N0 == "100.00":
                    burn(N0, S0, R0, loci)

    os.chdir(here)

