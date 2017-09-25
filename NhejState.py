import numpy as np
import os

class SimpleDsbState:
    def __init__(self):
        self.simpleDSB = 0
        self.ku_Li     = 0
        self.synapse   = 0
        self.null      = 0

    def initalize(self):
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
        if any(stateVector>1):
            raise ValueError("Invalid values for state")
    
