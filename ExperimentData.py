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

class Al_k_gamma:
	def __init__(self):
		self.data_27Gy  = np.array([[5,1,1.077],[7,0.782967,0.8516],[9.056,0.89011,0.9587],[11,0.725,0.8242],[13,0.6208,0.71978],[15,0.596154,0.642857],[17,0.516484,0.5879],[19,0.4423,0.4918],[21,0.346154,0.387363]])
		
		self.data_137Gy = np.array([[4.757,1.00,1.025],[5.2511,0.962887,0.9938],[7.15358,0.8619,0.9155],[9.239,0.8577,0.8907],[11.15,0.8268,0.8536],[13.148,0.7753,0.84536], 
						[15.220,0.6907,0.7485],[17.218,0.6639,0.7278],[19.132,0.6309,0.6948],[21.133,0.6227,0.6804],[23.127,0.5711,0.6227],[25.124,0.5320,0.6186], 
						[30.197,0.4330,0.4825],[35.110,0.3649,0.4371],[40.108,0.3134,0.3670],[45,0.2454,0.2969],[50,0.2268,0.2763],[55,0.2103,0.2536],[60,0.1938,0.2640]])

		self.data_27Gy[:,2] = self.data_27Gy[:,2] - self.data_27Gy[:,1]
		self.data_137Gy[:,2] = self.data_137Gy[:,2] - self.data_137Gy[:,1]
	
	def GetData_27Gy(self): 
		return self.data_27Gy
	
	def GetData_137Gy(self): 
		return self.data_137Gy