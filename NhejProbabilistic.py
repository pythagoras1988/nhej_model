import numpy as np 
import os 
import matplotlib.pyplot as plt
from NhejProbabilisticState import SimpleDsbState
from NhejProbabilisticState import ComplexDsbState
from NhejProbabilisticState import PairingStates

# *************************************************************************************************************************
# dsbMasterData format: N X 6 numpy Array 
# xPosition(Angstrom) X yPosition(Angstrom) X zPosition(Angstrom) X Complexity of break X Chromosome Index X Genomic Index
#**************************************************************************************************************************

global FLAG_configuration 
FLAG_configuration = False

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
		logdata = LogData()

		while(self.currTime < self.stopTime): 
			self.currTime += self.dt
			self._OneIteration()
			logdata.FillTime(self.currTime/60)
			logdata.CalculateRepairedState(self.stateList,'multiple')
			print('Current Time = %.1f mins...' %(self.currTime/60))

		logdata.ShowPlots()

	# Initialize DSBs in stateList. Each DSB gives rise to two states for the ends	
	def _InitializeDsbStates(self):
		count = 0
		for k in range(self.numDSB): 
			pos = self.data[k,0:2]
			gene_ID  = self.data[k,5] 
			chrom_ID = self.data[k,4]
			for kk in range(2):
				if self.data[k,3] == 1: 
					self.stateList.append(ComplexDsbState())
				else: 
					self.stateList.append(SimpleDsbState())
				self.stateList[-1].initialize(pos,gene_ID,chrom_ID)
				self.stateList[-1].ID = count # Set unique ID for each DSB
				count += 1
		self.numDSb = len(self.stateList) # update the size of numDSB to account for all ends

	# Initialize the Probability lists of all possible rejoining states
	def _InitializePairingStates(self): 	
		for k in range(self.numDSB-1): 
			for kk in range(k+1,self.numDSB):
				# The first Id s always smaller than the second ID
				self.pairingStateList.append(PairingStates())
				self.pairingStateList[-1].initialize(k,kk,self.stateList[k].startPosition,self.stateList[kk].startPosition)

		# Check the pairing state List contains the correct number of elements
		if len(self.pairingStateList)!=self.numPairingStates:
			raise Exceptions('Pairing State List construction errors...')

	def _OneIteration(self): 
		# Update Rejoined Probability for Pairing states
		for state in self.pairingStateList:
			state.RejoinedProb = self.ComputeRejoinProbability(state,'simple')

		# Compute Synapse formation probability for each DSB end
		synapseProbList = self.ComputeSynapseProbability(self.pairingStateList)

		# Update the state of each DSB end first
		counter = 0 
		for state in self.stateList: 
			state.stateUpdate(self.dt,synapseProbList[counter])
			counter += 1

	def ComputeSynapseProbability(self,pairingStateList): 
		# state is the pairingStateList
		probabilityList = []

		#use uniform distribution for prior!
		prior_pairedState = 1 / len(pairingStateList)

		for k in range(len(self.stateList)):
			tmpProb = 0
			for state in pairingStateList: 
				if state.ID_1==k or state.ID_2==k:
					tmpProb += state.RejoinedProb * prior_pairedState
			probabilityList.append(tmpProb)

		return probabilityList

	def ComputeRejoinProbability(self,state,option):
		# Use simple model for computing rejoin probability; Neglect boundary and centromere effect 
		if option == 'simple': 
			spatial_prob = 1 / np.sqrt(2*np.pi) * 0.5
			state_prob   = 1
			if NhejProcess.GetStateType(self.stateList[state.ID_1]) == 0: 
				state_prob *= self.stateList[state.ID_1].ku_XL 
			else:  
				state_prob *= self.stateList[state.ID_1].ku_PKcs_artemis

			if NhejProcess.GetStateType(self.stateList[state.ID_2]) == 0: 
				state_prob *= self.stateList[state.ID_2].ku_XL 
			else:  
				state_prob *= self.stateList[state.ID_2].ku_PKcs_artemis

			return spatial_prob*state_prob

	@staticmethod
	def GetStateType(state):
		stateName = state.__class__.__name__
		if stateName=='SimpleDsbState':
			return 0
		elif stateName=='ComplexDsbState':
			return 1
		else:
			raise Exception('Unknown States...')

class LogData:
	def __init__(self):
		self.dirName = 'Nhej_Repair_Outfiles'
		self.time = []
		self.plotData = []
		self.multiplePlotData = [[] for x in range(4)]

		if not os.path.isdir(self.dirName):
			os.mkdir(self.dirName)
	
	def FillTime(self,time): 
		self.time.append(time)

	def CalculateRepairedState(self,stateList,option): 
		if option == 'single': 
			self.plotData.append(stateList[0].null)
		elif option == 'mean': 
			tmp_sum = 0
			for k in range(len(stateList)): 
				tmp_sum += stateList[k].null 
			self.plotData.append(tmp_sum / (k+1) )
		elif option == 'multiple': 
			# only use this option if there are small number of DSB states
			for k in range(len(self.multiplePlotData)): 
				self.multiplePlotData[k].append(stateList[k].synapse)

	def ShowPlots(self): 
		for k in range(len(self.multiplePlotData)): 
			plt.plot(self.time,self.multiplePlotData[k],c=np.random.rand(3,1))
		#plt.plot(self.time,self.plotData,'r-')
		plt.xlabel('Time/mins')
		plt.ylabel('Probability')
		plt.show()





