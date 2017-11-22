import numpy as np
import os

class SimpleDsbState:
    def __init__(self):
        self.rateConstant = [1,1,1] #in per second
        self.simpleDSB = 0
        self.ku_Li     = 0
        self.synapse   = 0
        self.null      = 0
        self.ID.       = -1
        self.position  = np.array([0,0,0])
        self._initialize()

    def _initalize(self):
        self.simpleDSB += 1

    def getStateAsVector(self):
        return [self.simpleDSB,self.ku_Li,self.synapse,self.null]

    def stateChange1(self):
        self.simpleDSB -= 1
        self.ku_Li += 1

    def stateChange2(self):
        self.ku_Li -= 1
        self.synapse += 1

    def stateChange3(self):
        self.synpase -= 1
        self.null += 1

    def stateCheck(self):
        stateVector = self.getStateAsVector()
        if any(stateVector>1) or sum(stateVector)>1:
            raise ValueError("Invalid values for state")

    def stateStepping(self,dt):
        stateVector = self.getStateAsVector()
        ind = stateVector.index(1)
        if ind == 0 and np.random.uniform()<dt*self.rateConstant[0]:
            self.stateChange1
        elif ind == 1 and np.random.uniform()<dt*self.rateConstant[1]:
            self.stateChange2
        elif ind == 2 and np.random.uniform()<dt*self.rateConstant[2]:
            self.stateChange3
        else:
            pass

    def getDeltaTime(self):
        self.stateCheck()
        stateVector = self.getStateAsVector()
        ind = stateVector.index(1)
        r = np.random.uniform()
        return 1/self.rateConstant[ind] * np.log(1/r)

    def Set_ID(self,ID): 
        self.ID = ID

class ComplexUnjoinedDsbState:
    def __init__(self):
        self.rateConstant = [1,1]
        self.complexDSB = 0
        self.ku        = 0
        self.PKcs      = 0
        self.ID.       = -1
        self.position  = np.array([0,0,0])
        self._initialize()

    def _initialize(self):
        self.complexDSB += 1

    def getStateAsVector(self):
        return [self.complexDSB,self.ku,self.PKcs]

    def stateChange1(self):
        self.complexDSB -= 1
        self.ku += 1

    def stateChange2(self):
        self.ku -= 1
        self.PKcs += 1

    def stateCheck(self):
        stateVector = self.getStateAsVector()
        if any(stateVector>1) or sum(stateVector)>1:
            raise ValueError("Invalid values for state")

    def stateStepping(self,dt):
        stateVector = self.getStateAsVector()
        ind = stateVector.index(1)
        if ind == 0 and np.random.uniform()<dt*self.rateConstant[0]:
            self.stateChange1
        elif ind == 1 and np.random.uniform()<dt*self.rateConstant[1]:
            self.stateChange2
        else:
            pass

    def getDeltaTime(self):
        self.stateCheck()
        stateVector = self.getStateAsVector()
        ind = stateVector.index(1)
        r = np.random.uniform()
        return 1/self.rateConstant[ind] * np.log(1/r)

    def Set_ID(self,ID): 
        self.ID = ID

class ComplexJoinedDsbState:
    def __init__(self):
        self.rateConstant = [1,1,1]
        self.synapse = 0
        self.repairing = 0
        self.Li_XR = 0
        self.null = 0
        self.ID.  = -1
        self.position  = np.array([0,0,0])
        self._initialize()

    def _initialize(self):
        self.synapse += 1

    def getStateAsVector(self):
        return [self.synapse,self.reparing,self.Li_XR,self.null]

    def stateChange1(self):
        self.synapse -= 1
        self.repairing += 1

    def stateChange2(self):
        self.repairing -= 1
        self.Li_XR += 1

    def stateChange3(self):
        self.Li_XR -= 1
        self.null += 1

    def stateCheck(self):
        stateVector = self.getStateAsVector()
        if any(stateVector>1) or sum(stateVector)>1:
            raise ValueError("Invalid values for state")

    def stateStepping(self,dt):
        stateVector = self.getStateAsVector()
        ind = stateVector.index(1)
        if ind == 0 and np.random.uniform()<dt*self.rateConstant[0]:
            self.stateChange1
        elif ind == 1 and np.random.uniform()<dt*self.rateConstant[1]:
            self.stateChange2
        elif ind == 2 and np.random.uniform()<dt*self.rateConstant[2]:
            self.stateChange3
        else:
            pass

    def getDeltaTime(self):
        self.stateCheck()
        stateVector = self.getStateAsVector()
        ind = stateVector.index(1)
        r = np.random.uniform()
        return 1/self.rateConstant[ind] * np.log(1/r)

    def Set_ID(self,ID): 
        self.ID = ID
