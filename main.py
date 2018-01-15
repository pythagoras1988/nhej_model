##------------------------------------------------------------------------------------------------------
# This is the main file to be run to initialize NHEJ DNA repair module.
# The structure of the codes is as follows:
#								   _______________________________
# damageMat 1              -----> |	         NHEJ module 		  |
# 	   .				   -----> | a) classify DSB complexity    | ----> Predict Chromosome Aberration
#      . 				   -----> | b) system biology parameters  |	----> Predict Repair Rates
#      . 				   -----> | c) dsb ends diffusion settings| ----> Predict gamma H2AX formation
# damageMat 23             -----> | d) stopping condition         |
# chromosome position data -----> |_______________________________|
#
# Author: Higgsino
# Date: 23/8/2017
#--------------------------------------------------------------------------------------------------------


import numpy as np
import time
import os
from process_damage import ProcessDamage
from classify_damage import ClassifyDamage
from NhejProbabilistic import NhejProcess

global debug
debug = False

def Energy2Dose(energy): 
	mass = 1.0 * (1.5**3) * 10**-12 * 10**-3
	return energy* 1.6 * 10**-19 /(46*mass) 

class PlaceChromosome: 
	def __init__(self): 
		try: 
			self.chromPosArray = np.loadtxt('chromosome_positions_final.txt')
		except: 
			raise IOError('No Chromosome positions file present!!')

	def GetPosition(self): 
		return self.chromPosArray


if __name__=='__main__':
	#**********************************************************************
	# Input the damage files
	#**********************************************************************
	if debug:
		fname = 'damageMat0.txt'
		numDamageData = 1
		dsbMasterData = np.empty([0,6])
	else:
		## full cell nucleus irradiation; store all damageMat in 1 folder
		totalDose = 0
		dsbMasterData = np.empty([0,6])
		path = 'damageData'
		damageDirs = []

		## Get the number of damage data directories; There might be many photon damage data
		for root, dirs, files in os.walk(os.getcwd()+'/damageData'):
			for name in dirs: 
				print os.path.join(root,name)
				damageDirs.append(name)

	for dirName in damageDirs:	
		if len(damageDirs)>0:
			path = 'damageData/' + dirName
			dirs = os.listdir(path)
		else: 
			dirs = os.listdir(path)

		numDamageData = len(dirs) - 2  #1 of the file is edepMaster

		try:
			eDep = np.loadtxt(path+'/edepMaster.txt')
		except:
			raise IOError('No energy deposition files!!')

		#**********************************************************************
		# Run Algorithm to process the raw damageMat files and classify
		# the damage
		#
		# Raw damageMat ----> Process damage class ----> Classify damage class
		# 					(Enforce direct/indirect 	(det. ssb, dsb, direct
		#					   effect condition)		  indirect effect)
		#**********************************************************************

		for k in range(numDamageData):
			if not debug:
				totalDose += Energy2Dose(eDep[k])
				fname = path + '/damageMat' + str(k) + '.txt'
			procDamage = ProcessDamage(fname)

			if not procDamage.get_FLAG_NULLDATA():
				#if debug:
				#	np.savetxt('processed_damage_data.txt',procDamage.get_final_damage())
				classifyDamage = ClassifyDamage(procDamage.get_final_damage(),procDamage.get_base_damage())
				if debug:
					print classifyDamage.getDBSCAN_label()
					print classifyDamage.getDSB()
					print classifyDamage.getSSB()
					print classifyDamage.getDirectBreak()
					print classifyDamage.getIndirectBreak()
				if classifyDamage.getDSB()>=0:
					dsbMasterData = np.append(dsbMasterData,classifyDamage.getDSB_data(np.rint(22*np.random.uniform())),axis=0)
			print('Processing damage data number = %d; Dose = %.5f Gy...' %(k,totalDose))

	##************************************************************************************************************************
	# Run NHEJ repair code for DSB master data
	# 
	# dsbMasterData format: N X 6 numpy Array
	# xPosition(Angstrom) X yPosition(Angstrom) X zPosition(Angstrom) X Complexity of break X Chromosome Index x Genomic Index
	#*************************************************************************************************************************
		
	#-------------------------------------------------------------
	# Print Overall Statistics! 
	#-------------------------------------------------------------
	numSimpleBreaks  = sum(dsbMasterData[:,3]==0)
	numComplexBreaks = sum(dsbMasterData[:,3]==1)
	print('Total Number of DSBs = %d' %(len(dsbMasterData[:,0])))
	print('Simple DSB = %d, Complex DSB = %d' %(numSimpleBreaks,numComplexBreaks))
	time.sleep(3)
	nhejPtr = NhejProcess(dsbMasterData,PlaceChromosome().GetPosition())
