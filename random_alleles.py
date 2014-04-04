#!/usr/bin/env python

import numpy as np


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
