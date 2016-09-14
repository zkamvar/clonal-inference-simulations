#!/usr/bin/env python3.4

import numpy as np
import math as math
import sys

class snps:
	def __init__(self, nloc = 100, mu = 1e-5, nchrom = 10)
	self.nloc = nloc
	self.mu = mu
	self.nchrom = nchrom

class zk_locus:
	"""
	A Class representing a locus.
	holds:
		alleles
		frequencies
		mutation rate
	By default, this will generate a locus with six to ten alleles with a
	mutation rate of 1e-5.

		import random_alleles as ra
		locus1 = ra.zk_locus()
		locus2 = ra.zk_locus(nall = 23, mu = 1e-3)

		locus1.get_frequencies()
		locus2.get_frequencies()
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
		
		import random_alleles as ra
		nloc = 10
		alleles = [ra.zk_locus() for i in range(nloc)]
		loci = ra.zk_loci(alleles)
		loci.nloc()
		loci.get_locus_names()
		loci.get_frequencies()

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
