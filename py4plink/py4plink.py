import os
import sys
import numpy as np
import pandas as pd
import bitarray as ba
import gc

class Plink(object):


	def __init__(self, path=None,name=None, force=False):

		self.path=path
		self.name=name
		self.ext={'.fam','.bed','.bim'}
		self.bim=None
		self.bed=None
		self.fam=None
		self.N_probes=0
		self.n_probes=0
		self.N_ind=0
		self._currentSNP = None
		self.chunk_size=25000

		self._bedcode = {
			2: ba.bitarray('11'),
			9: ba.bitarray('10'), #TODO (high) NA data handle
			1: ba.bitarray('01'),
			0: ba.bitarray('00')
		}

		if force:
			self.read_fam()
			self.read_bed()
			self.read_bim()

	def read_fam(self):
		"""Read the FAM file to get information about individuals
		Family ID
		Individual ID
		Paternal ID
		Maternal ID
		Sex (1=male; 2=female; other=unknown)
		Label - phenotype should be in separate file
		"""

		self.fam = pd.read_table(os.path.join(self.path,self.name+'.fam'), sep=' ',
								 names= ['family', 'individual', 'paternal', 'maternal', 'sex', 'label'],
								 dtype={'names': ['family', 'individual', 'paternal', 'maternal', 'sex', 'label'],
										'formats': ['S', 'S', int, int, int, int]}
								 )
		self.N_ind=self.fam.shape[0]
		print('Number of Individuals: %d' % self.N_ind)


	def read_bim(self):

		self.bim = pd.read_table(os.path.join(self.path,self.name +'.bim'), sep='\t', header=None, names=['CHR', 'ID', 'distance', 'bp', 'allele1', 'allele2'],
								 dtype={'names': ['CHR', 'ID', 'distance', 'bp', 'allele1', 'allele2'],
										'formats': [int, 'S', int, int, 'S', 'S']}, iterator=True)

		self.N_probes = self.bim.shape[0]
		print('Number of Probes {} in {}'.format(self.N_probes, self.name + '.bim'))

	def get_fam(self):
		return self.fam

	def get_bim(self):
		return self.bim

	def get_bed(self):
		if self._currentSNP is None:
			self.read_bed()
		d=self.nextSNPs(self.N_probes)
		return d

	def read_bed(self):

		self._currentSNP=0
		self.bed = open(os.path.join(self.path,self.name+'.bed'), 'rb')
		magicNumber = ba.bitarray(endian="little")
		magicNumber.fromfile(self.bed, 2)
		mode = ba.bitarray(endian="little")
		mode.fromfile(self.bed, 1)
		K = (4 - self.N_ind % 4) if self.N_ind % 4 != 0 else 0
		self.nru = self.N_ind + K
		# check magic number
		if magicNumber != ba.bitarray('0011011011011000'):
			raise IOError("Magic number from Plink .bed file not recognized")

		if mode != ba.bitarray('10000000'):
			raise IOError("Plink .bed file must be in default SNP-major mode")


	def nextSNPs(self, b):
		'''
		Unpacks the binary array of genotypes and returns an n x b matrix of floats of
		normalized genotypes for the next b SNPs, where n := number of samples.

		Parameters
		----------
		b : int
			Number of SNPs to return.
		Returns
		-------
		X : np.array with dtype float64 with shape (n, b), where n := number of samples

		'''
		if self._currentSNP==self.N_probes:
			return None
		if self._currentSNP + b > self.N_probes:
			b=(self.N_probes-self._currentSNP)

		c = self._currentSNP
		n = self.N_ind
		slice= ba.bitarray(endian="little")

		bit_number=((2*(c+b)*self.nru)-(2*c*self.nru))/8
		slice.fromfile(self.bed, bit_number )

		X = np.array(slice.decode(self._bedcode), dtype="float64").reshape((b, self.nru)).T
		X = X[0:n, :]

		self._currentSNP += b
		gc.collect()
		return X.T
