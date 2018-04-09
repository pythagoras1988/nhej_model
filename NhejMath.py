import numpy as np
import scipy.special as sp

class Calculate_spatial_prob:
	# Set Static variables
	delta = 200 # in angstrom; distance for synapse formation
	D     = 5000 #100*100 # in angstrom^2/s
	def __init__(self,state,time,chromPosArr,option='full'):
		self.pos1 = state.pos1
		self.pos2 = state.pos2
		self.chromID1 = state.chrom_ID1
		self.chromID2 = state.chrom_ID2
		#self.D1       = 100*100 # in angstrom^2/s
		#self.D2       = self.D1/10
		self.time = time
		self.chromPosArr = chromPosArr
		self.prob = 0
		self.lowBound = 0
		self.upBound  = 0

		if option == 'debug':
			self.prob = self._CalculateProb_Debug()
		elif option == 'simple':
			self.prob = self._CalculateProb_Simple()
		else:
			self.prob = self._CalculateProb()

	def GetProbability(self):
		return self.prob

	def _CalculateProb_Simple(self):
		# determine global position
		self.pos1 = self.pos1 - 7500
		self.pos2 = self.pos2 - 7500
		minDist1  = 7500 - abs(self.pos1)
		minDist2  = 7500 - abs(self.pos2)
		#determine the global position
		self.pos1 = self.pos1 + self.chromPosArr[self.chromID1,:].T
		self.pos2 = self.pos2 + self.chromPosArr[self.chromID2,:].T

		if np.linalg.norm(self.pos1-self.pos2) < 3000:
			prob = 1.
			for k in range(3):
				#--------------------------------------------------------------------------------------------------------------
				# Model for scaling the diffusion constant to take into account the effect of chromosome territory and boundary
				#--------------------------------------------------------------------------------------------------------------
				sigma1 = np.sqrt(2*self.D*self.time)
				sigma2 = sigma1
				if sigma1>minDist1[k]:
					sigma1 = minDist1[k]
				if sigma2>minDist2[k]:
					sigma2 = minDist2[k]
				alpha1 = 1/(2*sigma1**2)
				beta1  = 1/(np.sqrt(2)*sigma2)
				#--------------------------------------------------------------------------------------------------------------
				gamma1 = (self.pos1[k]-self.pos2[k]+self.delta)/(np.sqrt(2)*sigma2)
				gamma2 = (self.pos1[k]-self.pos2[k]-self.delta)/(np.sqrt(2)*sigma2)
				prob *= 1./2* (sp.erf(gamma1*np.sqrt(alpha1/(alpha1+beta1**2)))-sp.erf(gamma2*np.sqrt(alpha1/(alpha1+beta1**2))))
			#print prob
		else:
			prob = 0.
		return prob


	def _CalculateProb_Debug(self):
		centre = np.array([0,0,0]) # for 1 chromosome scenario
		self.pos1 = self.pos1 - 7500
		self.pos2 = self.pos2 - 7500
		halfLength = 0.8 #in microns
		halfLength *= 10**4 # in angstrom
		prob = 1

		# Calculate the probability of intersection for x,y and z
		for k in range(3):
			self.lowBound = centre[k] - halfLength
			self.upBound  = centre[k] + halfLength
			prob *= self._ProdGaussian(self.pos1[k],self.pos2[k]) + self._ProdGaussian(self.pos1[k],2*halfLength - self.pos2[k]) \
			+ self._ProdGaussian(self.pos1[k],-2*halfLength - self.pos2[k]) + self._ProdGaussian(2*halfLength - self.pos1[k],self.pos2[k]) \
			+ self._ProdGaussian(2*halfLength - self.pos1[k], 2*halfLength - self.pos2[k]) + self._ProdGaussian(2*halfLength - self.pos1[k],-2*halfLength - self.pos2[k]) \
			+ self._ProdGaussian(-2*halfLength - self.pos1[k],self.pos2[k]) + self._ProdGaussian(-2*halfLength - self.pos1[k],2*halfLength - self.pos2[k]) \
			+ self._ProdGaussian(-2*halfLength - self.pos1[k],-2*halfLength - self.pos2[k])
			prob *= 2 * 250 * (4/3*np.pi)**(1/3)
		return prob

	def _CalculateProb(self):
		pass

	#calculate the integral of product of gaussian based on the mean position, B1 and B2
	def _ProdGaussian(self,B1,B2):
		# Assume similar standard deviation for both ends
		sig_1 = np.sqrt(2*self.D*self.time)
		sig_f = 1/np.sqrt(2)*sig_1
		mu_f  = (B1+B2)/2
		# Rescale the upper and lower bound
		upBound = (self.upBound - mu_f)/(np.sqrt(2)*sig_f)
		lowBound = (self.lowBound - mu_f)/(np.sqrt(2)*sig_f)

		return sig_f / (2*np.sqrt(2*np.pi)*sig_1**2) * (sp.erf(upBound)-sp.erf(lowBound))
