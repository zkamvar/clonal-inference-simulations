#!/usr/bin/env python3.4

import numpy as np
import math as math
import sys


class zk_locus:
	"""
	A Class representing a locus.
	holds:
		alleles
		frequencies
		mutation rate
	By default, this will generate a locus with six to ten alleles with a
	mutation rate of 1e-5.
	"""
	def __init__(self, nall = None, mu = None, amin = 6, amax = 10):
		if nall is not None:
			self.nall = nall
		else:
			self.nall = np.random.random_integers(amin, amax, 1)[0]
		if mu is not None:
			self.mu = mu
		else:
			self.mu = 1e-5
		self.alleles = self.new_alleles()
		self.freq = self.new_freqs()
		self.replen = int(self.get_alleles()[0][0])

	def get_frequencies(self):
		return(self.freq)

	def get_alleles(self):
		return(self.alleles)

	def get_mu(self):
		return(self.mu)

	def set_mu(self, mu):
		self.mu = mu

	def get_replen(self):
		return(self.replen)

	def new_alleles(self):
		nall = self.nall
		replen = np.random.random_integers(2, 6, 1)[0]
		allele_list = list()
		for i in range(nall):
			allele = str((replen * 100) + replen * (i + 1))
			allele_list.append(allele)
		return(allele_list)

	def new_freqs(self):
		nall = self.nall
		alleles = np.random.uniform(size = nall)
		allsum = sum(alleles)
		for i in range(0, len(alleles)):
			alleles[i] = alleles[i]/allsum
		return(alleles)


class zk_loci:
	"""
	A Class representing a set of loci.
	holds:
		a list of zk_locus objects
		mutation rates
		locus names
	
	This is a big wrapper for the zk_locus class, but this can return sensible
	information about your alleles. You can initialize it like so for 10 loci:
	
		nloc = 10
		alleles = [zk.locus() for i in range(nloc)]
		loci = zk_loci(alleles)

	And boom, there ya go
	"""
	def __init__(self, alleles):

		self.alleles = alleles
		self.locus_names = self.new_locus_names()

	def get_locus_names(self):
		return(self.locus_names)

	def nloc(self):
		return(len(self.locus_names))

	def nall(self):
		nalls = self.get_allele_names()
		nalls = [len(nalls[x]) for x in range(self.nloc())]
		return(nalls)

	def get_locus(self, index):
		return(self.alleles[index])
	
	def get_alleles(self):
		return(self.alleles)

	def get_allele_names(self):
		nloc = self.nloc()
		all_list = [self.get_locus(i).get_alleles() for i in range(nloc)]
		return(all_list)

	def get_allele_dict(self):
		alleles = self.get_alleles()
		lnames = self.get_locus_names()
		adict = dict()
		for i in range(len(lnames)):
			adict[lnames[i]] = alleles[i]
		return(adict)

	def get_frequencies(self):
		nloc = self.nloc()
		all_list = [self.get_locus(i).get_frequencies() for i in range(nloc)]
		return(all_list)		

	def get_mu(self):
		mu = [self.get_locus(i).get_mu() for i in range(self.nloc())]
		return(mu)

	def set_mu(self, index, mu):
		self.get_locus(index).set_mu(mu)

	def get_mu_dict(self):
		mu = self.get_mu()
		lnames = self.get_locus_names()
		mudict = dict()
		for i in range(len(lnames)):
			mudict[lnames[i]] = mu[i]
		return(mudict)

	def new_locus_names(self):
		locus_list = list()
		alist = self.alleles
		for i in range(len(alist)):
			replen = alist[i].get_replen()
			locus_specifier = str(i+1) + ":" + str(replen)
			locus_list.append("Locus " + locus_specifier)
		return locus_list

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
	locnames: an optional list of names for the loci.

Output: A histogram for each allele frequency.
'''
def plot_allele_probabilities(allele_probs, nall, locnames = None):
	namecount = 0
	for i in allele_probs:
		count = 0
		print("\\")
		if locnames is not None:
			print(" |" + locnames[namecount])
		else:
			print(" |")
		namecount = namecount + 1
		print("/")
		for j in range(0, nall):
			count = count + 1
			if nall > 9 and count < 10:
				strcount = " " + str(count)
			else:
				strcount = str(count)

			print(strcount+": "+"="*int(round(i[j]*100))+" "+str(round(i[j], 3)))
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

'''
Generate names for a specified number of loci.

Inputs:
	nloc: an integer specifying the number of loci.
	alist: a list of lists of allele names for each locus generated by 
		   get_allele_names. This is optional. 

Output:
	a list of locus names.
'''
def get_loci_names(nloc = 5, alist = None):
	locus_list = list()
	for i in range(nloc):
		if alist is not None:
			replen = alist[i][0][0]
			locus_specifier = str(i+1) + ":" + replen
		else:
			locus_specifier = str(i+1)
		locus_list.append("Locus " + locus_specifier)
	return locus_list

'''
Scale mutation rate by repeat length.

Since it's assumed that shorter repeat lengths mutate faster than longer ones,
we are going to scale a base mutation rate by raising it to the log of the 
repeat length:

mu^log(replen)

Inputs:
	mu: a mutation rate
	alist: a list of lists of allele names for each locus generated by
		   get_allele_names. This is required.
'''
def scale_mutation_rate(mu = 1e-5, alist = None):
	mout = list()
	for i in range(len(alist)):
		if alist is not None:
			replen = float(alist[i][0][0])
			pow = math.log(replen)
			mout.append(math.pow(mu, pow))
	return mout