# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd

class Data:
    def __init__(self, fileName):
        self.file = open(fileName)
        self.file.readline()
        self.file.readline()
        buff = self.file.readline().split()
        self.tN1 = buff[2]
        self.tN2 = buff[5]
        buff = self.file.readline().split()
        self.nm1 = int(buff[2])
        self.nm2 = int(buff[5])
        self.frame = int(self.file.readline().split()[2])
        self.file.readline()
    
    def getData(self):
        buff = self.file.readline()
        self.currTime = float(buff.split()[3])
        return pd.read_csv(self.file, nrows=self.nm1, header=None, delim_whitespace=True)
        
        
data = Data('Nr_SOL_Li+.dat')

#for i in range(data.frame):
a = data.getData()
b = data.getData()
