import numpy as np 
import os 
from NhejProbabilisticState import SimpleDsbState
from NhejProbabilisticState import ComplexDsbState

# *********************************************************************************************************
# dsbMasterData format: N X 5 Array
# xPosition(Angstrom) X yPosition(Angstrom) X zPosition(Angstrom) X Complexity of break X Chromosome Index
#**********************************************************************************************************

class NhejProcess: 
	def __init__(self,dsbMasterData): 