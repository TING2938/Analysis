import numpy as np
from glob import glob
import pandas as pd
import os

res = pd.DataFrame(columns=('ppm', 'Li/SOL', 'case', 'SOL_00', 'SOL_01', 'SOL_10', 'SOL_11',
                            'Li_00', 'Li_01', 'Li_10', 'Li_11'))

def getExcel(directory):    
    os.chdir(directory)
    noLi = True
    pwd = os.getcwd().split('\\')[-3:]
    for i in 'half', 'equal', 'twice', '4v1':
        if i in pwd:
            noLi = False
            break
    
    fname = ['Nr_SOL.dat', 'Nr_Li+.dat']
    
    
    r = []
    data =np.loadtxt(fname[0])
    data /= 0.01 * data.sum()
    r.append(data)
    data = np.loadtxt(fname[1])
    data /= 0.01 * data.sum()
    r.append(data)
    
    pwd = directory.split('\\')
    global res
    if not noLi:
        res = res.append({'ppm':pwd[-3], 'Li/SOL':pwd[-2], 'case':pwd[-1],\
                          'SOL_00':r[0][0], 'SOL_01':r[0][1], 'SOL_10':r[0][2], 'SOL_11':r[0][3],\
                          'Li_00':r[1][0], 'Li_01':r[1][1], 'Li_10':r[1][2], 'Li_11':r[1][3]},\
                          ignore_index=True)
    else:
        res = res.append({'ppm':pwd[-2], 'Li/SOL':"zero", 'case':pwd[-1],\
                                 'SOL-SOL':r[0], 'SOL-Li':np.nan, 'Li-SOL':np.nan}, ignore_index=True)
                  
def direct(rootDir):
    for root, dirs, files in os.walk(rootDir):
        for diri in dirs:
            if diri in "abc":
                getExcel(os.path.join(root, diri))

direct(r"D:\tmp\BmimTFSI_sol_Litfsi\BmimTFSI_water_LiTFSI\333K")

