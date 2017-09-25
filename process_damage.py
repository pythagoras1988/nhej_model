# Process the direct and indirect damage from damageMat files
# output a matrix indicating information of strands that are damaged

import numpy as np
import os.path

# This class will process the DNA damage from the raw output damage.mat file
class ProcessDamage:
    eThresh = 16.5 # in eV
    distThresh = 3.2 # in angstrom

    def __init__(self,fname):
        self.fname_read = fname
        self.data = np.array([])
        self.dataLen = 0
        self.FLAG_NULLDATA = False

        self.__readFile() # Read the data files into self.data
        if self.data.size != 0 :
            self.__get_total_damage()
        else:
            self.FLAG_NULLDATA = True

    def __readFile(self):
        try:
            self.data = np.loadtxt(self.fname_read)
            self.dataLen = len(self.data)
        except:
            raise Exception('Cannot load damageMat files!!')

    # Get final damage data after processing direct and indirect damage
    def get_final_damage(self):
        return self.data

    def get_base_damage(self):
        tmp = (self.data[:,1]==0)
        return self.data[tmp,:]

    def get_FLAG_NULLDATA(self):
        return self.FLAG_NULLDATA

    # main method to determine the final damage to the DNA; use as in input to the classification of damage step
    def __get_total_damage(self):

        # step 0: Remove Base dmage
        tmp = (self.data[:,1]==1)
        self.data = self.data[tmp,:]

        # step 1: remove damage which happens more than distance threshold for direct effect.
        #  This damage will not be considered at all!
        tmp  = (self.data[:,7]>self.distThresh) * (self.data[:,6]==0)
        tmp1 = (1-tmp).astype(bool)
        self.data = self.data[tmp1,:]
        self.dataLen = len(self.data[:,1])

        # step 2: Consolidate damage from direct damage especially for full edep data files.
        # This requires the direct damage data to be together, before sorting
        print self.data.shape
        self.data = self.__sum_direct_damage(self.data,self.dataLen)
        print self.data.shape

        # step 3: Remove direct damage withe eDep less than energy threshold
        tmp  = (self.data[:,8]<self.eThresh) * (self.data[:,6]==0)
        tmp1 = (1-tmp).astype(bool)
        self.data = self.data[tmp1,:]
        self.dataLen = len(self.data[:,0])

        # step 4: sort the data in increasing genomic index
        self.data = self.__insertionSort(self.data)
        np.savetxt('processed_damage_data.txt',self.data)

        # set NULLDATA FLAG if data is empty after post processing
        if self.data.size is 0:
            self.FLAG_NULLDATA = True

    # Same as consolidate damage in matlab version
    def __sum_direct_damage(self,data,length):
        if data.size == 0:
            return data
        else:
            tmp = data[1,:]
            tmpMat = np.zeros(length)
            tmpMat = tmpMat.astype(bool)
            tmpInd = 0
            for k in xrange(1,length):
                if data[k,6]==0 and all(data[k,0:3]==tmp[0:3]):
                    data[tmpInd,8] = data[tmpInd,8] + data[k,8]
                    tmpMat[k] = True
                    print('Combining direct damage data...')
                else:
                    tmp = data[k,:]
                    tmpInd = k
            tmp1 = (1-tmpMat).astype(bool)
            return data[tmp1,:]

    def __insertionSort(self,data):
        for k in xrange(1,self.dataLen):
            temp = data[k,:].copy()
            count = k
            while count>0 and temp[0]<data[count-1,0]:
                data[count,:] = data[count-1,:]
                count -= 1
            data[count,:] = temp
        return data
