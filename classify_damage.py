import numpy as np
from sklearn.cluster import DBSCAN

# Classify the damaged strand data into SSB, DSB and direct and indirect effect
class ClassifyDamage:
    def __init__(self,data,data_baseDamage):
	self.eps = 10
    	self.min_samples = 2
    	self.data = data
        self.data_GI = data[:,0] # only interested in the Genomic Index data
        self.data_GI_all = np.append(self.data_GI,data_baseDamage[:,0])
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
            data[k,3] = self._classify_dsb_complexity(k,tmp)
            data[k,4] = chromosome_index
        return data

    def _classify_dsb_complexity(self,k,tmp):
        threshold_distance = 10 #in term of bp
        # classify the complexity of DNA DSB. 1 represent complex and 0 represent Simple
        if sum(tmp)>2:
            # if more than 2 strand breaks in a DSB automatically labelled as complex
            return 1
        else:
            GI = self.data_GI[tmp]
            minGI = min(GI)
            maxGI = max(GI)
            minGI -= threshold_distance
            maxGI += threshold_distance
            counter = [1 for k in self.data_GI_all if (k<maxGI and k>minGI)]
            if sum(counter)>2:
                return 1
            else:
                return 0
