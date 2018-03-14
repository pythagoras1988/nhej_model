import numpy as np
from NhejProbabilistic import NhejProcess
from NhejProbabilisticState import SimpleDsbState
from NhejProbabilisticState import ComplexDsbState
from NhejProbabilisticState import PairingStates
from NhejMath import Calculate_spatial_prob

class OptimizerWrapper:
	def __init__(self,*argTuple):
		self.numTrials = 4
		#self.ParamList = np.array([10000])
		#self.ParamList = np.linspace(1000,10000,self.numTrials)
		#self.ParamList = np.linspace(0.02,0.1,5)/60 # For XL data
		#self.ParamList  = np.linspace(200,300,3)
		self.ParamList = np.array([4.5,30.0,60.,90.])/60

		for k in range(self.numTrials):
			fname = 'Nhej_Repair_Outfiles/DiffusionConstant' + str(k) + '.txt'
			#self._ChangeDiffusionConstant(k)
			self._ChangeRateConstSimple(k)
			self._ChangeRateConstComplex(k)
			#self._ChangeDeltaConstant(k)
			nhejPtr = NhejProcess(*argTuple)
			nhejPtr.SaveOptimData(fname)

	def _ChangeDiffusionConstant(self,k):
		Calculate_spatial_prob.D = self.ParamList[k]

	def _ChangeDeltaConstant(self,k):
		Calculate_spatial_prob.delta = self.ParamList[k]

	def _ChangeRateConstSimple(self,k):
		rateConstant = [self.ParamList[k],2.5155/60, 0.02/60]
		SimpleDsbState.rateConstant = rateConstant

	def _ChangeRateConstComplex(self,k):
		rateConstant = [self.ParamList[k],4.2257/60,4.5/60,2.7559/60, 0.02/60]
		ComplexDsbState.rateConstant = rateConstant
