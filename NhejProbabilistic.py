import numpy as np 
import os 
import math
import time
import matplotlib.pyplot as plt
import scipy.special as sp
from multiprocessing import Pool
from NhejProbabilisticState import SimpleDsbState
from NhejProbabilisticState import ComplexDsbState
from NhejProbabilisticState import PairingStates
from NhejMath import Calculate_spatial_prob
from ExperimentData import Cobalt_60_gamma
from ExperimentData import Al_k_gamma
from ExperimentData import Ce_137_gamma

# *************************************************************************************************************************
# dsbMasterData format: N X 6 numpy Array 
# xPosition(Angstrom) X yPosition(Angstrom) X zPosition(Angstrom) X Complexity of break X Chromosome Index X Genomic Index
#**************************************************************************************************************************

global FLAG_configuration 
global FLAG_saveData
global numProcess 
global numDSB_global 	
FLAG_configuration = False
FLAG_saveData = False
numProcess = 6
numDSB_global = 0

###---------------------------------------------------------------------------------------------------------------------------
# This portion of the code is for parallel Processing 																		 |
###---------------------------------------------------------------------------------------------------------------------------
def multi_process_wrapper(args): 																							#| 	
	return Calculate_spatial_prob(*args)																					#|	
																														    #| 	
def CalculateRejoinProbPerState(state,time,chromPosArr,stateList):															#|	
	delta = 200 # in angstrom; distance for synapse formation																#|	
	D = 100*100																												#|
    # determine global startPosition 																						#|
	pos1 = state.pos1 - 7500																								#|
	pos2 = state.pos2 - 7500																								#|
	minDist1  = 7500 - abs(pos1)																							#|
	minDist2  = 7500 - abs(pos2)																							#|
	#determine the global startPosition 																					#|
	pos1 = pos1 + chromPosArr[state.chromID1,:].T 																			#|
	pos2 = spos2 + chromPosArr[state.chromID2,:].T 																			#|
																															#|
	if np.linalg.norm(pos1-pos2) < 3000:																					#|
		prob = 1.																											#|
		for k in range(3):																									#|
			#-------------------------------------------------------------------------------------------------------------- #|
			# Model for scaling the diffusion constant to take into account the effect of chromosome territory and boundary #|
			#-------------------------------------------------------------------------------------------------------------- #|
			sigma1 = np.sqrt(2*D*time)
			sigma2 = sigma1
			if sigma1>minDist1[k]: 																							#|	
				sigma1 = minDist1[k]																						#|
			if sigma2>minDist2[k]:																							#|
				sigma2 = minDist2[k]																						#|
			alpha1 = 1/(2*sigma1**2)																						#|
			beta1  = 1/(np.sqrt(2)*sigma2)																					#|
			#-------------------------------------------------------------------------------------------------------------- #|
			gamma1 = (pos1[k]-pos2[k]+delta)/(np.sqrt(2)*sigma2)															#|
			gamma2 = (pos1[k]-pos2[k]-delta)/(np.sqrt(2)*sigma2)															#|
			prob *= 1./2* (sp.erf(gamma1*np.sqrt(alpha1/(alpha1+beta1**2)))-sp.erf(gamma2*np.sqrt(alpha1/(alpha1+beta1**2))))#|
		#print prob 																										#|
	else: 																													#|
		prob = 0. 																											#|
																															#|	
	stateProb = 1.																											#|
	try:																													#|
		state_prob *= self.stateList[state.ID_1].ku_XL 																		#|
	except:  																												#|	
		state_prob *= self.stateList[state.ID_1].ku_PKcs_artemis															#|
																															#|	
	try:																													#|
		state_prob *= self.stateList[state.ID_2].ku_XL 																		#|
	except:																													#|
		state_prob *= self.stateList[state.ID_2].ku_PKcs_artemis															#|
																															#|
	state.RejoinedProb = prob * stateProb 	

def ComputeSynapsePartialProb(stateList):
	partialArr = np.zeros(numDSB_global)
	for k in range(numDSB_global):
		for state in stateList: 
			if state.ID_1==k or state.ID_2==k : 
				partialArr[k] += state.RejoinedProb	
	return partialArr																										#|
##---------------------------------------------------------------------------------------------------------------------------

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
		self.stopTime = 10./60 # in hours
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
		#self._remove_small_fragments()
		logdata = LogData()

		##----------------------------------------
		# Start Multiprocessing 
		##----------------------------------------
		global numDSB_global
		numDSB_global = self.numDSB
		self.pool = Pool(processes=numProcess)

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
		print('Total Number of DSBs = %d ...' %self.numDSB)

	# Initialize the Probability lists of all possible rejoining states
	def _InitializePairingStates(self): 	
		for k in range(self.numDSB-1): 
			for kk in range(k+1,self.numDSB):
				# The first Id s always smaller than the second ID
				if self._RelevantPairingState(k,kk):
					self.pairingStateList.append(PairingStates())
					self.pairingStateList[-1].initialize(k,kk,self.stateList[k].startPosition,self.stateList[kk].startPosition, 
						self.stateList[k].chromosome_ID,self.stateList[kk].chromosome_ID)

		# Check the pairing state List contains the correct number of elements
		#if len(self.pairingStateList)!=self.numPairingStates:
		#	raise Exception('Pairing State List construction errors...')

		print('Number of Relevant Pairing States = %d ...' %len(self.pairingStateList))
		self.numPairingStates = len(self.pairingStateList)

	def _OneIteration(self): 
		# Update Rejoined Probability for Pairing states
		# use multiprocessing to speed up the process
		startTime = time.time()
		for state in self.pairingStateList:
			state.RejoinedProb = self.ComputeRejoinProbability(state)

		#for k in range(int(math.floor(self.numPairingStates/numProcess))):
		#	listOfTuples = [(self.pairingStateList[k],self.currTime,self.chromPos,self.stateList) for kk in range(k,k+numProcess)]
		#	self.pool.map(multi_process_wrapper,listOfTuples)
		print time.time() - startTime
		# Compute Synapse formation probability for each DSB end
		synapseProbList = self.ComputeSynapseProbability(self.pairingStateList)
		print time.time() - startTime

		# Update the state of each DSB end first
		counter = 0 
		for state in self.stateList: 
			state.stateUpdate(self.dt,synapseProbList[counter])
			counter += 1
			#print sum(state.getStateAsVector())

	def ComputeSynapseProbability(self,pairingStateList): 
		# state is the pairingStateList
		probabilityList = np.zeros(self.numDSB)

		#Divide pairingStateList for Multiprocessing
		parallelList = [ [] for k in range(numProcess)]
		counter = 0
		for state in pairingStateList:
			parallelList[counter%numProcess].append(state)
			counter += 1
		
		'''
		for k in range(len(self.stateList)):
			tmpProb = 0.
			#for state in pairingStateList: 
			#	if state.ID_1==k or state.ID_2==k:
			#		tmpProb += state.RejoinedProb		
			if (tmpProb>1.): 
				raise Exception('Pure Pairwise interaction hypothesis fails...')
			probabilityList[k] = tmpProb
		'''  
		return sum(self.pool.map(ComputeSynapsePartialProb,parallelList))

	# Calculate rejoining probability between any 2 dsb ends!
	# Rejoining Prob = spatial_intersection_prob x prob(end_1=ready) x prob(end_2=ready)
	def ComputeRejoinProbability(self,state):
		# Use simple model for computing rejoin probability; Neglect boundary and centromere effect 
		if self.rejoinModelOption == 'simple': 
			spatial_prob = Calculate_spatial_prob(state,self.D1,self.currTime,self.chromPos,'simple').GetProbability()
		# take into account boundary effect of the chromosome (single);debug purpose
		if self.rejoinModelOption == 'debug': 
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

	# Check for relevant pairing states; states that are far away will never rejoin together
	def _RelevantPairingState(self,k,kk):
		distanceThreshold = 3. # in Microns 
		distanceThreshold *= 10**4
		pos1 = self.stateList[k].startPosition - 7500. + self.chromPos[self.stateList[k].chromosome_ID,:].T
		pos2 = self.stateList[kk].startPosition - 7500. + self.chromPos[self.stateList[kk].chromosome_ID,:].T
		if np.linalg.norm(pos1-pos2)<distanceThreshold: 
			return True
		else: 
			return False

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
				if NhejProcess.GetStateType(stateList[k]): 
					#complex
					tmp_sum += stateList[k].null
					#tmp_sum += stateList[k].ku_PKcs_artemis
					#tmp_sum += stateList[k].ku_PKcs
					#tmp_sum += stateList[k].ku_PKcs_XL
					#tmp_sum += stateList[k].synapse	
				else: 
					#simple
					tmp_sum += stateList[k].null
					#tmp_sum += stateList[k].ku_XL
					#tmp_sum += stateList[k].synapse
			self.plotData.append(tmp_sum/(k+1.))
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
		data = Al_k_gamma().GetData_27Gy()
		data = Ce_137_gamma().GetData_40Gy()
		#plt.errorbar(data[:,0],data[:,1],yerr=data[:,2],label='Al_k')
		plt.plot(data[:,0]*60,data[:,1]/max(data[:,1]),'ro',label='Ce-137')

class ChromosomeAberrationAlgo: 
	def __init__(self):
		pass
		







