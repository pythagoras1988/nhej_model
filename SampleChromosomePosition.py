import numpy as np 
import os

## Create an array of chromosome positions in the nucleus. The nucleus is centered at [0,0,0]

class SampleChromosomePosition:
	def __init__(self): 
		bufferDistance = 1000 # in angstrom
		self.numChromosome = 46
		self.radiusNucleus = 6 # in microns 
		self.radiusNucleus *= 10**4 # convert to angstroms
		self.lengthChromosome = 1.5 # in microns
		self.lengthChromosome *= 10**4 # convert to angstroms
		self.overlapDistance = self.lengthChromosome-bufferDistance
		self.chromPosArr = np.zeros((self.numChromosome,3))
		self._GenerateRandomChromosome()
		self._SaveAsAscii()

	def GetPositionArrays(self): 
		return self.chromPosArr

	@staticmethod
	def ReadChromosomeAsciiFile(): 
		return np.loadtxt('Chromosome_positions_Arrays.txt')

	def _SaveAsAscii(self): 
		np.savetxt('Chromosome_positions_Arrays.txt',self.chromPosArr)

	def _GenerateRandomChromosome(self):
		maxLengthChromosome = np.sqrt(3*self.lengthChromosome/2)
		count = 0
		while self._CheckOverlap() and count<100000:
			for k in range(self.numChromosome):
				x = (2*np.random.uniform()-1)*(self.radiusNucleus-maxLengthChromosome)
				y = (2*np.random.uniform()-1)*(self.radiusNucleus-maxLengthChromosome)
				z = (2*np.random.uniform()-1)*(self.radiusNucleus-maxLengthChromosome)
				self.chromPosArr[k,:] = np.array([x,y,z])
			print('Trial number %d ...' %count)
			count += 1

	def _CheckOverlap(self):
		# overlap between chromosomes
		for k in range(self.numChromosome-1): 
			for kk in range(k+1,self.numChromosome):
				distX = abs(self.chromPosArr[k,0] - self.chromPosArr[kk,0])
				distY = abs(self.chromPosArr[k,1] - self.chromPosArr[kk,1])
				distZ = abs(self.chromPosArr[k,2] - self.chromPosArr[kk,2])

				if distX<self.overlapDistance and distY<self.overlapDistance and distZ<self.overlapDistance:
					return True
				
		return False

if __name__ == '__main__': 
	SampleChromosomePosition().GetPositionArrays()
