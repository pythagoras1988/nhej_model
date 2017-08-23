import numpy as np 
from sklearn.cluster import DBSCAN

# Classify the damaged strand data into SSB, DSB and direct and indirect effect
class ClassifyDamage: 
	def __init__(self,data): 
		self.eps = 10 
    	self.min_samples = 2
    	self.data = data 
    	assert self.data.size == 0

    	#implement clustering via DBSCAN
    	self.label = self.__startDBSCAN() 

    def __startDBSCAN(self): 
    	# return cluster labels for the final damage matrix 
    	clusterPtr = DBSCAN(self.eps,self.min_samples)
    	return clusterPtr.fit_predict(self.data)

    def getDSB(self): 
    	return max(self.label)

    def getSSB(self):
    	return sum(self.label==-1)

    def getDirectBreak(self):
    	tmp = (data[:,6]==0 * data[:,1]==1) 
    	return sum(tmp)

    def getIndirectBreak(self):
    	tmp = (data[:,6]==1 * data[:,1]==1) 
    	return sum(tmp)

    def getDSB_data(self): 
    	pass
