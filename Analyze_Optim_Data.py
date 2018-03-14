import numpy as np
from ExperimentData import Ce_137_gamma
import matplotlib.pyplot as plt

class AnalyzeData:
    def __init__(self):
        self.dirName ='Nhej_Repair_Outfiles/InitialBindingConstant'
        self.numFiles =1
        self._PlotSimulationData()
        self._PlotData()

    def _PlotSimulationData(self):
        for k in range(self.numFiles):
            fname = self.dirName + str(k) + '.txt'
            data = np.loadtxt(fname)
            plt.plot(data[0,:],data[1,:],c=np.random.rand(3,1),label = str(k))

    def _PlotData(self):
        data = Ce_137_gamma().GetData_40Gy()
        plt.plot(data[:,0]*60,data[:,1]/30.,'ro',label='Ce-137')
        plt.xlabel('Time/mins')
        plt.xlim(0,10)
        plt.ylabel('Probability')
        plt.legend()
        plt.show()

if __name__ == '__main__':
    AnalyzeData()
