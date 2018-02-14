import numpy as np
from NhejProbabilistic import NhejProcess
from NhejProbabilisticState import SimpleDsbState
from NhejProbabilisticState import ComplexDsbState
from NhejProbabilisticState import PairingStates
from NhejMath import Calculate_spatial_prob

class OptimizerWrapper:
	def __init__(self,*argTuple):
		self.numTrials = 5
		#self.ParamList = np.linspace(10000,20000,self.numTrials)
		self.ParamList = np.linspace(0.02,1.0,5)/60

		for k in range(self.numTrials):
			fname = 'Nhej_Repair_Outfiles/XLConstant' + str(k) + '.txt'
			#self._ChangeDiffusionConstant(k)
			self._ChangeRateConstSimple(k)
			self._ChangeRateConstComplex(k)
			nhejPtr = NhejProcess(*argTuple)
			nhejPtr.SaveOptimData(fname)

	def _ChangeDiffusionConstant(self,k):
		Calculate_spatial_prob.D = self.ParamList[k]

	def _ChangeDeltaConstant(self,k):
		Calculate_spatial_prob.delta = self.ParamList[k]

	def _ChangeRateConstSimple(self,k):
		rateConstant = [1.,2.5155/60, self.ParamList[k]]
		SimpleDsbState.rateConstant = rateConstant

	def _ChangeRateConstComplex(self,k):
		rateConstant = [1.,4.2257/60,4.5/60,2.7559/60, self.ParamList[k]]
		ComplexDsbState.rateConstant = rateConstant
