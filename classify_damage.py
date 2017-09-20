import numpy as np
from sklearn.cluster import DBSCAN

# Classify the damaged strand data into SSB, DSB and direct and indirect effect
class ClassifyDamage:
    def __init__(self,data):
	self.eps = 10
    	self.min_samples = 2
    	self.data = data
        self.data_GI = data[:,0] # only interested in the Genomic Index data
    	assert self.data.size != 0

    	#implement clustering via DBSCAN
    	self.label = self._startDBSCAN()

    def _startDBSCAN(self):
    	# return cluster labels for the final damage matrix
    	clusterPtr = DBSCAN(self.eps,self.min_samples).fit(self.data_GI.reshape(-1,1))
    	return clusterPtr.labels_

    def getDBSCAN_label(self):
        return self.label

    def getDSB(self):
    	return max(self.label)

    def getSSB(self):
    	return sum(self.label==-1)

    def getDirectBreak(self):
    	tmp = (self.data[:,6]==0) * (self.data[:,1]==1)
    	return sum(tmp)

    def getIndirectBreak(self):
    	tmp = (self.data[:,6]==1) * (self.data[:,1]==1)
    	return sum(tmp)

    #**************************************************************************
    # This method extracts the state of the DSB including:                    #
    # 1) Complexity of breaks                                                 #
    # 2) Position of breaks                                                   #  
    #                                                                         #  
    # Col 1,2,3: X,Y,Z Position
    # Col 4 : Simple or complex break (1 is complex)
    # Col 5 : Chromosome Index
    #**************************************************************************
    def getDSB_data(self,chromosome_index):
        #Initialization
    	data = np.zeros((self.getDSB(),5))
        for k in xrange(self.getDSB()+1):
            tmp    = (self.label==k)
            tmpDat = self.data[tmp,:]
            data[k,0] = np.mean(tmpDat[:,3]) # x position
            data[k,1] = np.mean(tmpDat[:,4]) # y position
            data[k,2] = np.mean(tmpDat[:,5]) # z position
            data[k,3] = self._classify_dsb_complexity()
            data[k,4] = chromosome_index

    def _classify_dsb_complexity(self):
        pass
