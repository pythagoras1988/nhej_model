import numpy as np
import os

class SimpleDsbState:
	#------------------------------------------------------------------------------------------------
	# Simple DSB end ---(rate_1)---> Binding of Ku and XL ---(diffusion)---> Synapsis Formation
	# ---(rate_2)---> ku release ---(rate_3)---> XL released (rejoined or null state)
	#------------------------------------------------------------------------------------------------
    # rate constant list as class static variable
    rateConstant = [1.,2.5155/60, 0.05/60] #in per second
    def __init__(self):
        self.simpleDSB = 0
        self.ku_XL     = 0
        self.synapse   = 0
        self.XL        = 0
        self.null      = 0
        self.genomic_ID     = -1
        self.chromosome_ID  = -1
        self.synapse_ID    = -1
        self.ID        = -1
        self.startPosition  = np.array([0,0,0])

    def initialize(self,pos,gene_ID,chrom_ID):
        self.simpleDSB += 1
        self.startPosition = pos
        self.genomic_ID = gene_ID
        self.chromosome_ID = chrom_ID

    def getStateAsVector(self):
        return [self.simpleDSB,self.ku_XL,self.synapse,self.XL,self.null]

    def formSynapse(self):
    	if self.synapse>0:
    		return true

    def stateUpdate(self,dt,synapseProb):
        tmpVector = self.getStateAsVector()
    	self.simpleDSB -= self.simpleDSB*self.rateConstant[0]*dt
    	self.ku_XL += tmpVector[0]*self.rateConstant[0]*dt - synapseProb*self.ku_XL
        self.synapse +=  synapseProb*tmpVector[1] - self.synapse*self.rateConstant[1]*dt
    	self.XL = self.XL + tmpVector[2]*self.rateConstant[1]*dt - self.XL*self.rateConstant[2]*dt
    	self.null += tmpVector[3]*self.rateConstant[2]*dt

    def stateCheck(self):
    	# Probability should sum up to 1
    	#if sum(self.getStateAsVector()!=1):
         #   raise ValueError("Sum of Probability not equal to 1")
        if any(i > 1 for i in self.getStateAsVector()) or any(i < 0 for i in self.getStateAsVector()):
        	raise ValueError("Invalid probability value!")


class ComplexDsbState:
	#------------------------------------------------------------------------------------------------
	# Complex DSB end ---(rate_1)---> Binding of Ku,DNA-Pkcs,Artemis ---(diffusion)---> Synapsis Formation
	# ---(rate_2)---> artemis processing ---(rate_3)---> ku_DNA-PKcs_XL ---(rate_4)---> XL ---(rate_5)---> (rejoined or null state)
	#------------------------------------------------------------------------------------------------
    # Rate constant as class static variable
    rateConstant = [1.,4.2257/60,4.5/60,2.7559/60, 0.05/60] # in seconds
    def __init__(self):
        self.complexDSB = 0
        self.ku_PKcs_artemis = 0
        self.synapse     = 0
        self.ku_PKcs     = 0
        self.ku_PKcs_XL  = 0
        self.XL          = 0
        self.null        = 0
        self.genomic_ID     = -1
        self.chromosome_ID  = -1
        self.synapse_ID     = -1
        self.ID          = 0
        self.startPosition  = np.array([0,0,0])

    def initialize(self,pos,gene_ID,chrom_ID):
        self.complexDSB += 1
        self.startPosition = pos
        self.genomic_ID = gene_ID
        self.chromosome_ID = chrom_ID

    def getStateAsVector(self):
        return [self.complexDSB,self.ku_PKcs_artemis,self.synapse,self.ku_PKcs,self.ku_PKcs_XL,self.XL,self.null]

    def formSynapse(self):
    	if self.synapse>0:
    		return true

    def stateUpdate(self,dt,synapseProb):
        tmpVector = self.getStateAsVector()
    	self.complexDSB -= tmpVector[0]*self.rateConstant[0]*dt
    	self.ku_PKcs_artemis += tmpVector[0]*self.rateConstant[0]*dt - synapseProb*tmpVector[1]
        self.synapse += synapseProb*tmpVector[1] - tmpVector[2]*self.rateConstant[1]*dt
    	self.ku_PKcs += tmpVector[2]*self.rateConstant[1]*dt - tmpVector[3]*self.rateConstant[2]*dt
        self.ku_PKcs_XL += tmpVector[3]*self.rateConstant[2]*dt - tmpVector[4]*self.rateConstant[3]*dt
        self.XL += tmpVector[4]*self.rateConstant[3]*dt - tmpVector[5]*self.rateConstant[4]*dt
    	self.null += tmpVector[5]*self.rateConstant[4]*dt

    def SetRateConstant(self,rateConstant):
        if not isinstance(rateConstant,list) or len(rateConstant)!=5:
            raise Exception('Wrong Inputs during optimization')
        else:
            self.rateConstant = rateConstant

    def stateCheck(self):
    	# Probability should sum up to 1
    	#if sum(self.getStateAsVector()!=1):
         #   raise ValueError("Sum of Probability not equal to 1")
        if any(i > 1 for i in self.getStateAsVector()) or any(i < 0 for i in self.getStateAsVector()):
        	raise ValueError("Invalid probability value!")

class PairingStates:
    def __init__(self):
        self.ID_1 = -1
        self.ID_2 = -1
        self.chrom_ID1 = -1
        self.chrom_ID2 = -1
        self.pos1 = np.array([0,0,0])
        self.pos2 = np.array([0,0,0])
        self.RejoinedProb = 0.

    def initialize(self,id1,id2,pos1,pos2,chrom_ID1,chrom_ID2):
        self.ID_1 = id1
        self.ID_2 = id2
        self.pos1 = pos1
        self.pos2 = pos2
        self.chrom_ID1 = chrom_ID1
        self.chrom_ID2 = chrom_ID2
