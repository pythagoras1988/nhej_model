import numpy as np 
import os 
from NhejProbabilisticState import SimpleDsbState
from NhejProbabilisticState import ComplexDsbState
from NhejProbabilisticState import PairingStates

# *************************************************************************************************************************
# dsbMasterData format: N X 6 numpy Array 
# xPosition(Angstrom) X yPosition(Angstrom) X zPosition(Angstrom) X Complexity of break X Chromosome Index X Genomic Index
#**************************************************************************************************************************

class NhejProcess: 
	def __init__(self,dsbMasterData): 
		self.currTime = 0 # in seconds
		self.stopTime = 3 # in hours
		self.stopTime *= 3600 # in seconds
		self.dt       = 0.5 # in seconds
		self.D1       = 170*100 # in angstrom^2/s
		self.D2       = self.D1/10
		self.data     = dsbMasterData 
		self.numDSB   = 0
		self.stateList = []
		self.pairingStateList = []

		self.numDSB   = len(self.data[:,0])
		self.numPairingStates = 0.5 * self.numDSB * (self.numDSB-1)
		self._InitializeDsbStates()
		self._InitializePairingStates()

		while(self.currTime < stopTime): 
			self._OneIteration()
			self.currTime += self.dt
			print('Current Time = %.1f mins...' %(self.currTime/60))

	# Initialize DSBs in simple and complex lists		
	def _InitializeDsbStates(self):
		for k in range(self.numDSB): 
			pos = self.data[k,0:2]
			gene_ID  = self.data[k,5] 
			chrom_ID = self.data[k,4]
			if self.data[k,3] == 1: 
				self.stateList.append(ComplexDsbState())
			else: 
				self.stateList.append(SimpleDsbState())
			stateList[-1].initialize(pos,gene_ID,chrom_ID)
			stateList[-1].ID = k # Set unique ID for each DSB

	# Initialize the Probability lists of all possible rejoining states
	def _InitializePairingStates(self): 	
		for k in range(self.numDSB-1): 
			for kk in range(k+1,self.numDSB):
				# The first Id s always smaller than the second ID
				self.pairingStateList.append(PairingStates())
				self.pairingStateList[-1].initialize(k,kk,self.stateList[k].startPosition,self.stateList[kk].startPosition)

		# Check the pairing state List
		if len(self.pairingStateList)!=self.numPairingStates:
			raise Exceptions('Pairing State List construction errors...')

	def _OneIteration(self): 






