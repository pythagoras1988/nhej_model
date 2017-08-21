import numpy as np 
from sklearn.cluster import DBSCAN

# Classify the damaged strand data into SSB, DSB and direct and indirect effect
class ClassifyDamage: 
	def __init__(self,data): 
		self.eps = 10 
    	self.min_samples = 2
    	self.data = data 
    	#implement clustering via DBSCAN
    	self.__startDBSCAN() 

    def __startDBSCAN(self): 
    	clusterPtr = DBSCAN(self.eps,self.min_samples)
    	return clusterPtr.fit_predict(self.data)

    def getDSB(self): 
    	pass
    def getSSB(self):
    	pass
    def getDirectBreak(self):
    	pass
    def getIndirectBreak(self):
    	pass
