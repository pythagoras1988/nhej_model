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

global debug 
debug = True

if __name__=='__main__':
	#**********************************************************************
	# Input the damage files 
	#**********************************************************************
	if debug: 
		fname = 'damageMat0.txt'
	else: 
		## To be fixed! 
		fname = 'damageMat0.txt'

	#**********************************************************************
	# Run Algorithm to process the raw damageMat files and classify 
	# the damage
	#
	# Raw damageMat ----> Process damage class ----> Classify damage class
	# 					(Enforce direct/indirect 	(det. ssb, dsb, direct
	#					   effect condition)		  indirect effect) 	
	#**********************************************************************

	procDamage = ProcessDamage(fname)
	if not procDamage.get_FLAG_NULLDATA(): 
		#if debug: 
		#	np.savetxt('processed_damage_data.txt',procDamage.get_final_damage())	
		classifyDamage = ClassifyDamage(procDamage.get_final_damage())



