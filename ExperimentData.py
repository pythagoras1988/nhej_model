import numpy as np 

# From Sternalow 2001
class Cobalt_60_gamma: 
	def __init__(self): 
		# Time in mins
		self.residual_dsb = np.array([[7.97849,0.903746,0.937587],[23.0323,0.587697,0.672298],[42,0.391878,0.539086],[65.9355,0.319699,0.383998]])
		self.residual_dsb[:,2] = self.residual_dsb[:,2]-self.residual_dsb[:,1]
		# Time in hours
		self.num_dsb = np.array([[0.131753,514.943,533.333],[0.380924,335.632,383.142],[0.693933,223.755,306.513],[1.09376,182.375,217.625],[3.05513,55.1724,76.6284],[5.88113,41.3793,52.1073],[22.5576,24.5211,42.9119]])
		self.num_dsb[:,2] = self.num_dsb[:,2] - self.num_dsb[:,1]
		self.num_dsb[:,[1,2]] = self.num_dsb[:,[1,2]] / 565

	def GetResidualDSB(self): 
		return self.residual_dsb

	def GetNumDSB(self): 
		return self.num_dsb

