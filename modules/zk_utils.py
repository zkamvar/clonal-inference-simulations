#!/usr/bin/env python3

import random_alleles as ra
import numpy as np
import pandas as pd


'''
Parser for input parameters for the argparse module.

Parameters:
    arg_line a string containing a parameter and its value[s]

Output:
    a generator containing the parameter and its value[s]

Example:
    
    res = convert_arg_line_to_args("--sexrate 0.0 1.0")
    [x for x in res]
    ## ['--sexrate', '0.0', '1.0']

    # It's designed to be used with argparse:

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        fromfile_prefix_chars = '@'
        )

    parser.convert_arg_line_to_args = convert_arg_line_to_args
'''
def convert_arg_line_to_args(arg_line):
    for arg in arg_line.split():
        if not arg.strip():
            continue
        if arg[0] == '#':
            break
        yield arg

'''
Write arguments to a file for use with argparse

Parameters:
    args a dict object of arguments produced from argparse. Note: this must
    contain an agrument called "cfg" specifying the name of the configuration
    file to write to.
Output:
    a new file with a different argument on each line

Example:
    
    # Assuming you have already set up your parser:

    pars = parser.parse_args()
    args_to_file(pars)

'''
def args_to_file(args):
    f = open(args.cfg, "w")
    for arg, value in vars(args).items():
        if not isinstance(value, (list)):
            value = [value]
        val = '{} '*len(value)
        f.write("--{} ".format(arg) + val.format(*value) + "\n")
    f.close

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

'''
A wrapper function for the random_alleles module that will generate a zk_loci
object

Parameters:
    nloc - the number of loci to include
    murate - the mutation rate (must be a list of length nloc)
    amax - the maximum number of alleles/locus
    amin - minimum number of alleles/locus

Output:
    A zk_loci object.

Examples:
    
    generate_loci(5, [5e-5]*5, 2, 20)
'''
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