import numpy as np
from NhejProbabilistic import NhejProcess
from NhejProbabilisticState import SimpleDsbState
from NhejProbabilisticState import ComplexDsbState
from NhejProbabilisticState import PairingStates
from NhejMath import Calculate_spatial_prob

class OptimizerWrapper:
	def __init__(self,*argTuple):
		self.numTrials = 3
		self.ParamList = np.linspace(10000,20000,self.numTrials)

		for k in range(self.numTrials):
			fname = 'Nhej_Repair_Outfiles/DiffusionConstant' + str(k) + '.txt'
			self._ChangeDiffusionConstant(k)
			nhejPtr = NhejProcess(*argTuple)
			nhejPtr.SaveOptimData(fname)

	def _ChangeDiffusionConstant(self,k):
		Calculate_spatial_prob.D = self.ParamList[k]

	def _ChangeDeltaConstant(self,k):
		Calculate_spatial_prob.delta = self.ParamList[k]

	def _ChangeRateConstSimple(self,k):
		rateConstant = []
		SimpleDsbState.rateConstant = rateConstant

	def _ChangeRateConstComplex(self,k):
		rateConstant = []
		ComplexDsbState.rateConstant = rateConstant
