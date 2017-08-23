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

        self.__readFile(self.fname_read) # Read the data files into self.data
        if self.data.size not 0: 
            self.__get_total_damage()
        else: 
            self.FLAG_NULLDATA = True

    def __readFile(self):
        try:
            os.pathisfile(self.fname_read):
            self.data = np.loadtxt(self.fname_read)
            self.dataLen = len(self.data)
        except:
            raise Exception('Cannot load damageMat files!!')

    # Get total energy deposited in 1 event
    def get_energy_deposit(self):
        pass

    # Get final damage data after processing direct and indirect damage
    def get_final_damage(self):
        return self.data

    def get_FLAG_NULLDATA(self): 
        return self.FLAG_NULLDATA 

    # main method to determine the final damage to the DNA; use as in input to the classification of damage step
    def __get_total_damage(self):
        # sort the data in increasing genomic index
        self.data = self.__insertionSort(self.data)

        # Consolidate damage from direct damage especially for full edep data files
        self.data = self.__sum_direct_damage(self.data,self.dataLen)

        # remove damage which happens more than distance threshold and less than energy threshold for direct damage
        tmp  = (self.data[:,7]>disThresh * self.data[:,6]==0 * self.data[:,8]<eThresh)
        tmp1 = tmp * np.arange(self.dataLen)
        tmp1 = tmp1[tmp.astype(bool)]
        self.data = np.delete(self.data,tmp1,axis=0)

        # set NULLDATA FLAG if data is empty after post processing
        if self.data.size is 0: 
            self.FLAG_NULLDATA = True

    # Same as consolidate damage in matlab version
    def __sum_direct_damage(self,data,length):
        if data.size ==0:
            return data
        else:
            tmp = data[1,:]
            tmpMat = np.zeros(length,1)
            tmpMat.astype(bool)
            tmpInd = 0

            for k in xrange(1,length):
                if data[k,6]==0 and (data[k,0:2]==tmp[1,0:2]).all()
                    data[tmpInd,8] = data[tmpInd,8] + data[k,8]
                    tmpMat[k] = True
                else:
                    tmp = data[k,:]
                    tmpInd = k
            tmp1 = tmpMat * np.arange(length)
            tmp1 = tmp1[tmpMat]
            return np.delete(data,tmp1,axis=0)

    def __insertionSort(self,data):
        for k in xrange(1,self.dataLen):
            tmp = data[k,:]
            count = k
            while count>0 and tmp[0]<data[count-1,1]:
                data[count,:] = data[count-1,:]
                count = count - 1
            data[count,:] = tmp

