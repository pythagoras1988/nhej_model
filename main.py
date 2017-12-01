##-------------------------------------------------------------------------
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
#--------------------------------------------------------------------------


import numpy as np
from process_damage import ProcessDamage
from classify_damage import ClassifyDamage
from NhejProcess import NhejPtr

global debug
debug = True

if __name__=='__main__':
	#**********************************************************************
	# Input the damage files
	#**********************************************************************
	if debug:
		fname = 'damageMat0.txt'
		numDamageData = 1
		dsbMasterData = np.empty([0,5])
	else:
		## full cell nucleus irradiation; store all damageMat in 1 folder
		totalEnergy = 0
		dsbMasterData = np.array([])
		path = '/damageData/'
		dirs = os.listdir(path)
		numDamageData = len(dirs) - 1 #1 of the file is edepMaster
		try:
			eDep = np.loadtxt(path+'edepMaster.txt')
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
			totalEnergy += eDep[k]
			fname = dirs(k)
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
				#dsbMasterData.append(classifyDamage.getDSB_data(np.rint(23*np.random.uniform())))
				dsbMasterData = np.append(dsbMasterData,classifyDamage.getDSB_data(np.rint(23*np.random.uniform())),axis=0)


	##**********************************************************************
	# Run NHEJ repair code for DSB master data
	#***********************************************************************
	nhejPtr = NhejPtr(dsbMasterData)
