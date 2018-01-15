import numpy as np 
import os 
import matplotlib.pyplot as plt
from NhejProbabilisticState import SimpleDsbState
from NhejProbabilisticState import ComplexDsbState
from NhejProbabilisticState import PairingStates
from NhejMath import Calculate_spatial_prob
from ExperimentData import Cobalt_60_gamma

# *************************************************************************************************************************
# dsbMasterData format: N X 6 numpy Array 
# xPosition(Angstrom) X yPosition(Angstrom) X zPosition(Angstrom) X Complexity of break X Chromosome Index X Genomic Index
#**************************************************************************************************************************

global FLAG_configuration 
global FLAG_saveData
FLAG_configuration = False
FLAG_saveData = False

class NhejProcess: 
	def __init__(self,dsbMasterData,chromPos): 
		##-----------------------------------------
		# Program I/O 
		##-----------------------------------------
		self.plotOption = 'mean'
		self.rejoinModelOption = 'simple'
		##-----------------------------------------
		# Program Model Parameters
		##-----------------------------------------		
		self.currTime = 0.1 # in seconds
		self.stopTime = 30./60 # in hours
		self.stopTime *= 3600 # in seconds
		self.dt       = 0.5 # in seconds
		self.D1       = 100*100 # in angstrom^2/s
		self.D2       = self.D1/10
		self.data     = dsbMasterData 
		self.chromPos = chromPos
		self.numDSB   = 0
		self.stateList = []
		self.pairingStateList = []

		self.numDSB   = len(self.data[:,0])
		self.numPairingStates = 0.5 * 2 * self.numDSB * (2*self.numDSB-1)
		self._InitializeDsbStates()
		self._InitializePairingStates()
		logdata = LogData()

		while(self.currTime < self.stopTime): 
			self.currTime += self.dt
			self._OneIteration()
			logdata.FillTime(self.currTime/60)
			logdata.CalculateRepairedState(self.stateList,self.plotOption,0)
			if FLAG_saveData: 
				timeInterval = 0 
				if self.currTime>timeInterval:
					logdata.SaveSynapseRepairData(self.pairingStateList,self.currTime)
					timeInterval += 0.5 # Save the rejoined data after every 0.5 seconds

			print('Current Time = %.2f mins...' %(self.currTime/60))

		logdata.ShowPlots()

	# Initialize DSBs in stateList. Each DSB gives rise to two states for the ends	
	def _InitializeDsbStates(self):
		count = 0
		for k in range(self.numDSB): 
			pos = self.data[k,0:3]
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
		self.numDSB = len(self.stateList) # update the size of numDSB to account for all ends

	# Initialize the Probability lists of all possible rejoining states
	def _InitializePairingStates(self): 	
		for k in range(self.numDSB-1): 
			for kk in range(k+1,self.numDSB):
				# The first Id s always smaller than the second ID
				self.pairingStateList.append(PairingStates())
				self.pairingStateList[-1].initialize(k,kk,self.stateList[k].startPosition,self.stateList[kk].startPosition, 
					self.stateList[k].chromosome_ID,self.stateList[kk].chromosome_ID)

		# Check the pairing state List contains the correct number of elements
		if len(self.pairingStateList)!=self.numPairingStates:
			raise Exception('Pairing State List construction errors...')

	def _OneIteration(self): 
		# Update Rejoined Probability for Pairing states
		for state in self.pairingStateList:
			state.RejoinedProb = self.ComputeRejoinProbability(state,self.rejoinModelOption)
			#print state.RejoinedProb

		# Compute Synapse formation probability for each DSB end
		synapseProbList = self.ComputeSynapseProbability(self.pairingStateList)

		# Update the state of each DSB end first
		counter = 0 
		for state in self.stateList: 
			state.stateUpdate(self.dt,synapseProbList[counter])
			counter += 1
			#print sum(state.getStateAsVector())

	def ComputeSynapseProbability(self,pairingStateList): 
		# state is the pairingStateList
		probabilityList = []

		#use uniform distribution for prior!
		prior_pairedState = 1. / self.numDSB

		for k in range(len(self.stateList)):
			tmpProb = 0.
			for state in pairingStateList: 
				if state.ID_1==k or state.ID_2==k:
					tmpProb += state.RejoinedProb
			if (tmpProb>1.): 
				raise Exception('Pure Pairwise interaction hypothesis fails...')
			probabilityList.append(tmpProb)

		return probabilityList

	# Calculate rejoining probability between any 2 dsb ends!
	# Rejoining Prob = spatial_intersection_prob x prob(end_1=ready) x prob(end_2=ready)
	def ComputeRejoinProbability(self,state,option):
		# Use simple model for computing rejoin probability; Neglect boundary and centromere effect 
		if option == 'simple': 
			spatial_prob = Calculate_spatial_prob(state,self.D1,self.currTime,self.chromPos,'simple').GetProbability()
		# take into account boundary effect of the chromosome (single);debug purpose
		if option == 'debug': 
			spatial_prob = Calculate_spatial_prob(state,self.D1,self.currTime,self.chromPos,'debug').GetProbability()
			
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

# Class to save data into ascii files and output to stdout
class LogData:
	def __init__(self):
		self.dirName = 'Nhej_Repair_Outfiles'
		self.option = None
		self.time = []
		self.plotData = []
		self.multiplePlotData = [[] for x in range(4)]
		self.stateList_simple = [[] for x in range(5)]
		self.stateList_complex = [[] for x in range(7)]

		if not os.path.isdir(self.dirName):
			os.mkdir(self.dirName)

		# Create files for saving data
		if FLAG_saveData:
			self.repairData_fname = self.dirName + '/' + 'Synapse_Repair_Data.txt'
			if os.path.isfile(self.repairData_fname):
				os.remove(self.repairData_fname)
	
	def FillTime(self,time): 
		self.time.append(time)

	def CalculateRepairedState(self,stateList,option,ID=0): 
		self.option = option
		if option == 'single': 
			if NhejProcess.GetStateType(stateList[ID]): 
				# Complex
				self.stateList_complex[0].append(stateList[ID].complexDSB)
				self.stateList_complex[1].append(stateList[ID].ku_PKcs_artemis)
				self.stateList_complex[2].append(stateList[ID].synapse)
				self.stateList_complex[3].append(stateList[ID].ku_PKcs)
				self.stateList_complex[4].append(stateList[ID].ku_PKcs_XL)
				self.stateList_complex[5].append(stateList[ID].XL)
				self.stateList_complex[6].append(stateList[ID].null)
			else:
				# Simple
				self.stateList_simple[0].append(stateList[ID].simpleDSB)
				self.stateList_simple[1].append(stateList[ID].ku_XL)
				self.stateList_simple[2].append(stateList[ID].synapse)
				self.stateList_simple[3].append(stateList[ID].XL)
				self.stateList_simple[4].append(stateList[ID].null)

		elif option == 'mean': 
			tmp_sum = 0
			for k in range(len(stateList)): 
				tmp_sum += stateList[k].null 
			self.plotData.append(1 - tmp_sum / (k+1) )
		elif option == 'multiple': 
			# only use this option if there are small number of DSB states
			for k in range(len(self.multiplePlotData)): 
				self.multiplePlotData[k].append(stateList[k].null)

	def ShowPlots(self): 
		#print self.multiplePlotData
		if self.option == 'multiple':
			for k in range(len(self.multiplePlotData)): 
				plt.plot(self.time,self.multiplePlotData[k],c=np.random.rand(3,1),label = str(k))
		if self.option == 'single':
			if not self.stateList_simple[0]: 
				for k in range(len(self.self.stateList_complex)): 
					plt.plot(self.time,self.stateList_complex[k],c=np.random.rand(3,1),label = str(k))
			else: 
				for k in range(len(self.stateList_simple)): 
					plt.plot(self.time,self.stateList_simple[k],c=np.random.rand(3,1),label = str(k))
		if self.option == 'mean':
			# Plot the remaining DSBs over time
			plt.plot(self.time,self.plotData,'r',label='DSBs Percentage')

		self._PlotExperimentData()
		plt.xlabel('Time/mins')
		plt.ylabel('Probability')
		#plt.yscale('log')
		plt.xlim(0,75)
		plt.legend()
		plt.show()

	def SaveSynapseRepairData(self,pairingStateList,time):
		# Only save non-Zero probabilities
		with open(self.repairData_fname,'a') as writeFile: 
			writeFile.write('%.1f \t' %time) # in seconds
			for state in pairingStateList: 
				if state.RejoinedProb > 0: 
					writeFile.write('%d \t %d \t %.5f \n' %(state.ID_1, state.ID_2, state.RejoinedProb)) 

	def _PlotExperimentData(self):
		data = Cobalt_60_gamma().GetNumDSB()
		plt.errorbar(data[:,0]*60,data[:,1],yerr=data[:,2],label='Co-60')

class ChromosomeAberrationCalc: 
	def __init__(self):
		pass
		







