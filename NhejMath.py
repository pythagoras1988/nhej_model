import numpy as np 
import scipy.special as sp

class Calculate_spatial_prob:
	def __init__(self,state,D,time,option='full'): 
		self.pos1 = state.pos1 
		self.pos2 = state.pos2 
		self.chromID1 = state.chrom_ID1 
		self.chromID2 = state.chrom_ID2
		self.D = D
		self.time = time
		self.prob = 0

		if option == 'debug': 
			self.prob = self._CalculateProbDebug
		else: 
			self.prob = self._CalculateProb

	def GetProbability(self): 
		return self.prob

	def _CalculateProbDebug(self): 
		centre = np.array([0,0,0]) # for 1 chromosome scenario
		halfLength = 0.8 #in microns
		halfLength *= 10**4 # in angstrom
		lowBound = centre - halfLength
		upBound  = centre + halfLength

	def _CalculateProb(self): 
		pass

	#calculate the integral of product of gaussian based on the mean position, B1 and B2
	def _ProdGaussian(B1,B2,upBound,lowBound):
		# Assume similar standard deviation for both ends
		sig_1 = sqrt(2*self.D*self.time)
		sig_f = 1/sqrt(2)*sig_1
		mu_f  = (B1+B2)/2
		# Rescale the upper and lower bound
		upBound = (upBound - mu_f)/(sqrt(2)*sig_f)
		lowBound = (lowBound - mu_f)/(sqrt(2)*sig_f)

		return sig_f / (2*sqrt(2*np.pi)*sig_1**2)*(sp.erf(upBound)-sp.erf(lowBound))