#!/usr/bin/env python
def write(N0, C0, S0, R0):
    with open("./clone_"+str(N0)+"_pop_"+str(S0)+"_rep_%02d" % R0+".cfg", 'w') as g:
        g.write(
            "[optimized]\n\n"
            "[S0] # Population Size\n"
            "S0="+str(S0)+"\n\n"
            "[Sam0] # Sample Size\n"
            "Sam0=[10, 25, 50, 100]\n\n"
            "[G0] # Number of Generations\n" 
            "G0=10000\n\n"
            "[C0] # Percent Clonal Reproduction\n"
            "C0="+'%.2f' % C0+"\n\n"
            "[R0] # Number of Replicates\n"
            "R0=10\n\n" 
            "[L0] # Number of Loci\n"
            "L0=[1]*10\n\n"
            "[N0] # Fileset Name\n"
            "N0=clone_"+str(N0)+"_pop_"+str(S0)+"_rep_%02d" % R0+"\n"
            )
        g.close()


def burn(N0, S0, R0):
    with open("./BURNIN_pop_"+str(S0)+"_rep_%02d" % R0+".cfg", 'w') as g:
        g.write(
            "[optimized]\n\n"
            "[S0] # Population Size\n"
            "S0="+str(S0)+"\n\n"
            "[B0] # Number of Burnin Generations\n" 
            "B0=1000\n\n"
            "[L0] # Number of Loci\n"
            "L0=[1]*10\n\n"
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
#sizes = (10000, 10000)
clones = (100, 99.99, 99.95, 99.9, 99, 95, 90, 80, 50, 0)
#for S0 in sizes:
for C0 in clones:
    N0 = name(C0)
    print(N0+"\n")
    for R0 in range(100):

        '''
        if R0 < 10:
            print("%02d" % R0)
            s = "clone_"+str(N0)+"_pop_10000_rep_%02d" % R0+".cfg"
            print(s.split("_")[5].split(".")[0])
        '''

        write(N0, C0, 10000, R0)
        if N0 == "100.00":
            burn(N0, 10000, R0)
        
