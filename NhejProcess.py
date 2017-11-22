import numpy as np
import os
import NhejState


class NhejPtr:
	def __init__(self,data):
		self.totalTime = 4 * 3600 # in seconds
		self.currTime  = 0
		self.numDSB    = len(data[:,0])
		self.D1        = 170*100 # in angstrom^2/s
		self.D2        = self.D1/10
		self.data = data
		self.stateList = [] # list containing state instance objects
		self.repairedList = [] # list containing all states that have done repairing 

		# Initialize all dsb states
		self._initializeDsbStates()

		# Run each iterations using Gillespie algorithm
		if self.currTime<self.totalTime or numDSB==0:
			self.OneIteration()
			self.PrintInfo()
			# log data!

	##-------------------------------------------------------------------
	# Initialize initial state and position
	###------------------------------------------------------------------
	def _initializeDsbStates(self):
		count = 0;
		for k in range(self.data[:,0]):
			if dsbMasterData[k,3] == 1:
				self.stateList.append(NhejState.ComplexUnjoinedDsbState())
				self.stateList[-1].position = self.data[k,0:2]
				self.stateList[-1].Set_ID(count)	
				self.stateList.append(NhejState.ComplexUnjoinedDsbState())
				self.stateList[-1].position = self.data[k,0:2]
				self.stateList[-1].Set_ID(count)
				count += 1	
			else:
				self.stateList.append(NhejState.SimpleDsbState())
				self.stateList[-1].position = self.data[k,0:2]
				self.stateList[-1].Set_ID(count)
				count += 1	

		# Update total number of breaks
		self.numDSB = len(self.stateList)

	##-------------------------------------------------------------------
	# Update state and position in 1 time step 
	###------------------------------------------------------------------
	def OneIteration(self):
		dt = self._timeStepping()
		self.currTime += dt
		print("Current time is %f mins",%(self.currTime/60))
		for k in range(numDSB):
			#Update the state in 1 time step
			self.stateList[k].stateStepping()
			#Update the position in 1 time step
			if self.GetStateType(self.stateList[k])<2: 
				self.stateList[k].position = self.stateList[k].position + self.FindStepSize(self.D2,dt)
			else:
				self.stateList[k].position = self.stateList[k].position + self.FindStepSize(self.D1,dt)
		self._Process_Rejoining()

	##-------------------------------------------------------------------
	# Algorithm to check for DSB rejoining 
	###------------------------------------------------------------------	
	def _Process_Rejoining(self): 
		pass

	##-------------------------------------------------------------------
	# Basic class to support the iteration class above
	###------------------------------------------------------------------

	def _timeStepping(self):
		dt_List = [self.stateList[k].getDeltaTime() for k in self.numDSB]
		return min(dt_List)

	def FindStepSize(self,D,t):
		stepVector = np.array([0,0,0])
		std = np.sqrt(2*D*t)
		stepVector[0] = np.random.normal(0.0,std)
		stepVector[1] = np.random.normal(0.0,std)
		stepVector[2] = np.random.normal(0.0,std)
		return stepVector

	def GetStateType(self,state): 
		stateName = state.__class__.__name__
		if stateName=='SimpleDsbState': 
			return 0 
		elif stateName=='ComplexUnjoinedDsbState': 
			return 1 
		elif stateName=='ComplexJoinedDsbState': 
			return 2 
		else: 
			raise Exception('Unknown States...')

	##-------------------------------------------------------------------
	# Destructor
	###------------------------------------------------------------------
	def PrintInfo(self): 
		def __init__(self): 
			self._printStats()

		def _printStats(self):
			print('-------------------------------------------------')
			print('Current Time = %f mins....' %self.currTime/60)
			print('Total number of unprocessed DSBs = %d ....' %self.numDSB)
			print('-------------------------------------------------')


