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
		self.maxIndex  = 0 
		self.data = data
		self.stateList = [] # list containing state instance objects
		self.rejoinedData = [] # list containing information of the DSBs that have finished rejoining
		self.repairedList = [] # list containing all states that have done repairing 

		# Initialize all dsb states
		self._initializeDsbStates()

		# Run each iterations using Gillespie algorithm
		if self.currTime<self.totalTime or numDSB==0:
			self.OneIteration()
			self.PrintInfo()
		# log data!
		logdata = LogData()
		logdata.Save_Rejoined_Data(self.rejoined_data)
		logdata.Save_Repaired_List(self.repairedList)

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

		self.maxIndex = count
		# Update total number of breaks
		self.numDSB = len(self.stateList)

	##-------------------------------------------------------------------------------------------------------------------
	# Update state and position in 1 time step 
	# 
	# 		1. Determine time step : dt 
	#					|
	#					|
	#					V
	#		2. Determine state probability and transition : stateStepping method in NhejState Class
	#      				|
	#					|
	#					V
	#       3. Determine position of each state after diffusion in dt time: _timeStepping method and FindStepSize method
	# 					|
	#					|
	#					V
	#		4. Determine which 2 DSB state will rejoin by forming a synapse: _Process_Rejoin method
	#					|
	#					|
	# 					V
	# 		5. Determine which rejoined state will be fully repaired: _Process_Repair method
	###-------------------------------------------------------------------------------------------------------------------
	def OneIteration(self):
		minTimeStep = 0.5 # Minimum time step
		dt = self._timeStepping()
		if dt > minTimeStep:
			dt = minTimeStep # in seconds
		self.currTime += dt
		print("Current time is %f mins",%(self.currTime/60))
		for k in range(self.numDSB):
			#Update the state in 1 time step
			self.stateList[k].stateStepping()
			#Update the position in 1 time step
			if self.GetStateType(self.stateList[k])<2: 
				self.stateList[k].position = self.stateList[k].position + self.FindStepSize(self.D2,dt)
			else:
				self.stateList[k].position = self.stateList[k].position + self.FindStepSize(self.D1,dt)
		self._CheckHitBoundary()
		self._Process_Rejoin()
		self._Process_Repair()

	##-------------------------------------------------------------------
	# Algorithm to check for DSB rejoining 
	###------------------------------------------------------------------	
	import sklearn.metrics

	def _Process_Rejoin(self): 
		# Set Rejoining condition
		rejoin_distance = 10 # in angstrom
		rejoin_prob  = 0.5 # Probabilitistic rejoining

		# initialize position arrays
		positionArrays = np.zeros((self.numDSB,3))
		for k in range(self.numDSB): 
			positionArrays[k,:] = self.stateList[k].position

		# Compute pairwise distance matrix
		#pw_distance_matrix = self._Compute_Pairwise_Distance(positionArrays)

		# Calculate which pair of DSBs will rejoin
		rejoined_list = [] 
		for k in range(0,self.numDSB-1): 
			for kk in range(k+1,self.numDSB): 
				dist = np.linalg.norm(positionArrays[k,:]-positionArrays[kk,:]) 
				rand = np.random.uniform()
				if (dist<=rejoin_distance and rand<=rejoin_prob): 
					rejoined_list.append(k) 
					rejoined_list.append(kk) 
					print('|-------------------------------------------|')
					print('|DSB rejoining happens between %d and %d ...|' %(k,kk))
					print('|-------------------------------------------|')

					# push rejoined state to the stateList
					tmpPosition = (stateList[k].position() + stateList[kk].position())/2 
					stateList.append(NhejState.ComplexJoinedDsbState()) 
					stateList[-1].position = tmpPosition 
					stateList[-1].set_ID(self.maxIndex)
					self.maxIndex += 1

					# Store rejoined state information into rejoinedList
					tmp_rejoinedState = [None] * 8 # initialization
					tmp_rejoinedState[0] = self.currTime 
					tmp_rejoinedState[1] = tmpPosition[0]
					tmp_rejoinedState[2] = tmpPosition[1]
					tmp_rejoinedState[3] = tmpPosition[2]
					tmp_rejoinedState[4] = k 
					tmp_rejoinedState[5] = self.GetStateType(self.stateList[k])
					tmp_rejoinedState[6] = kk 
					tmp_rejoinedState[7] = self.GetStateType(self.stateList[kk])
					self.rejoinedData.append(tmp_rejoinedState) #Append the information of rejoining into rejoinedData!!

		# delete rejoined DSBs from the stateList 
		del stateList[rejoined_list]
		self.numDSB = len(self.stateList) #Update number of states in the list

	def _Compute_Pairwise_Distance(self,positionArray): 
		pw_distance_matrix = sklearn.metrics.pairwise.pairwise_distances(positionArray,metric='euclidean')
		return pw_distance_matrix

	##-------------------------------------------------------------------
	# Class to process the completion of repair of DSB 
	###------------------------------------------------------------------
	def _Process_Repair(self): 
		tmp_repairedList = []
		for k in range(len(self.numDSB)): 
			# Check that the DSB has finished repairing
			if self.GetStateType(self.stateList[k])==2: 
				if self.stateList[k].null==1: 
					# Push the repaired state information to the repairedList
					tmp_repairedState = [None] * 4 # initialization
					tmp_repairedState[0] = self.currTime 
					tmp_repairedState[1] = self.stateList[k].position[0]
					tmp_repairedState[2] = self.stateList[k].position[0]
					tmp_repairedState[3] = self.stateList[k].position[0]

					tmp_repairedList.append(k)
					self.repairedList.append(tmp_repairedState) 
		del self.stateList[tmp_repairedList]
		self.numDSB = len(self.stateList)

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
	# Method to handle exception when position stepping reaches the  
	# boundary of the cell nucleus 
	###------------------------------------------------------------------
	def _CheckHitBoundary(self): 
		pass

	##-------------------------------------------------------------------
	# class to save data to ascii file 
	###------------------------------------------------------------------

	class LogData: 
		def __init__(self): 
			self.dirName = 'Nhej_Repair_Outfiles'
			if not os.path.isdir(self.dirName): 
				os.mkdir(self.dirName)

		def Save_Rejoined_Data(self,data): 
			# Data is a list data with N x 8 elements
			with open(self.dirName + '/rejoined_data.txt') as saveFile: 
				# Write header
				saveFile.write('Time/mins \t xPos/A \t yPos/A \t zPos/A \t index_1 \t Type_1 \t index_2 \t Type_2 \n')
				for k in range(len(data)): 
					saveFile.write('%.2f \t %.1f \t %.1f \t %.1f \t %d \t %d \t %d \t %d \n' 
						%(data[k][0],data[k][1],data[k][2],data[k][3],data[k][4],data[k][5],data[k][6],data[k][7]))

		def Save_Repaired_List(self,data): 
			# Data is a list data with N x 4 elements
			with open(self.dirName + '/repaired_data.txt') as saveFile: 
				# Write header
				saveFile.write('Time/mins \t xPos/A \t yPos/A \t zPos/A \n')
				for k in range(len(data)): 
					saveFile.write('%.2f \t %.1f \t %.1f \t %.1f \n' 
						%(data[k][0],data[k][1],data[k][2],data[k][3]))
	##-------------------------------------------------------------------
	# Destructor
	###------------------------------------------------------------------
	class PrintInfo: 
		def __init__(self): 
			self._printStats()

		def _printStats(self):
			print('-------------------------------------------------')
			print('Current Time = %f mins....' %self.currTime/60)
			print('Total number of unprocessed DSBs = %d ....' %self.numDSB)
			print('-------------------------------------------------')


