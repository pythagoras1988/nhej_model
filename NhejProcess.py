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
		self.stateList = [] # contain state objects

		# Initialize all dsb states
		self._initializeDsbStates()

		# Run each iterations using Gillespie algorithm
		if self.currTime<self.totalTime:
			self.oneIteration()
			# log data!

	def _initializeDsbStates(self):
		for k in range(data[:,0]):
			if dsbMasterData[k,3] == 1:
				self.stateList.append(NhejState.ComplexUnjoinedDsbState())
			else:
				self.stateList.append(NhejState.SimpleDsbState())

	def oneIteration(self):
		dt = self._timeStepping()
		self.currTime += dt
		print("Current time is %f mins",%(self.currTime/60))
		for k in range(numDSB):
			self.stateList[k].stateStepping()

	def _timeStepping(self):
		dt_List = [self.stateList[k].getDeltaTime() for k in self.numDSB]
		return min(dt_List)

	def positionStepping(self):
		pass
