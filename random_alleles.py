#!/usr/bin/env python

import numpy as np
import sys


def rescale_allele_probabilities(allele_probs, nloc, nall):
	for i in range(0, nloc):
		allsum = sum(allele_probs[i])
		for j in range(0, nall):
			allele_probs[i][j] = allele_probs[i][j]/allsum
	return allele_probs

'''
Generate allele probabilities from a random uniform distribution for a certain
number of alleles. 

Inputs:
	nloc: an integer. The number of loci.
	nall: an integer. The number of alleles per locus. 
Output:
	a numpy array with the number of rows equal to the number of loci.
	each row sums to one.
'''
def get_allele_probabilities(nloc, nall):
	allele_probs = np.random.uniform(size = nloc*nall).reshape(nloc, nall)
	allele_probs = rescale_allele_probabilities(allele_probs, nloc, nall)
	return allele_probs

'''
Plot the allele probabilities on the command line.

Inputs:
	allele_probs: a numpy array where each row sums to 1.
	nall: the number of columns. 

Output: A histogram for each allele frequency.
'''
def plot_allele_probabilities(allele_probs, nall):
	for i in allele_probs:
		count = 0
		print("|")
		for j in range(0, nall):
			count = count + 1
			print(str(count)+": "+"="*int(round(i[j]*100))+" "+str(round(i[j], 3)))
	return

'''
Generate names of microsatellite alleles.

Inputs:
	nloc: an integer defining the number of loci.
	nall: an integer defining the number of alleles per locus.

Output:
	A list of lists of allele names. Each locus will have a separate repeat
	length defined randomly. This repreat length will be reflected in the allele
	names. So, if the repeat is of length 3 (CAT)^n, then the allele names will
	be '3XX' and they will be divisible by 3.
'''
def get_allele_names(nloc, nall):
	replens = np.random.random_integers(2, 6, nloc)
	locus_list = list()
	for i in range(len(replens)):
		prefix = replens[i]
		allele_list = list()
		for j in range(nall):
			allele = str((prefix * 100) + prefix * (j + 1))
			sys.stdout.write(allele + " ")
			allele_list.append(allele)
		locus_list.append(allele_list)
		print("\n")
	return locus_list
