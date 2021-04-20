import numpy as np
from glob import glob
import matplotlib.pyplot as plt
import os

def plotN(directory):
    global index1, index2, fn, title, fig, ax, ddt, tt1, tt2, ll
    
    os.chdir(directory)
    noLi = True
    pwd = directory.split('\\')[-3:]
    for i in 'half', 'equal', 'twice', '4v1':
        if i in pwd:
            noLi = False
            break
    
    if noLi:
        ndx1 = 0
    else:
        ndx1 = index1[pwd[-2]]
    
    try:
        fileName = glob(fn)[0]
        data = np.loadtxt(fileName)
        x = data[1:, 0]
        y = data[1:, 1] * 100 / sum(data[1:, 1])
        
        t1 = int(x[0])
        t2 = int(x[-1])
        l = t2 - t1
        if l > 30: dt = 10
        elif l > 10: dt = 5
        elif l >= 5: dt = 2
        else: dt = 1 
        
        if ll <= l:
            ll = l
            ddt = dt
            tt1 = t1
            tt2 = t2
        
        ax[ndx1].plot(x, y, "-o", color='#FF6400', lw=1.5)
        ax[ndx1].bar(x, y, color='#00A383', width = 0.6)
        #ax[ndx1][ndx2].set_title(title[j])
        
        ax[ndx1].set_ylabel('(%)')
    except:
        ax[ndx1].axis('off')


fig, ax = plt.subplots(5, 1,  figsize=(3, 9), sharex=True)

allfn = 'Energy-2-1*.dat', 'Energy-2-3*.dat', 'Energy-4-1*.dat', 'Energy-4-3*.dat', \
            'Energy-4-4*.dat', 'Energy-4-2*.dat', 'Energy-2-4*.dat'
alltitle = 'Li-Bmim', 'Li-Tf2N', 'SOL-Bmim', 'SOL-Tf2N', 'SOL-SOL', 'SOL-Li', 'Li-SOL'

ff = 4

fn = allfn[ff]
title = alltitle[ff]
ddt = 0
tt1 = 0
tt2 = 0
ll = 0

ax[4].set_xlabel('N (#)')
index1 = {'zero':0, 'half':1, 'equal':2, 'twice':3, '4v1':4}

def directory(rootDir):
    for root, dirs, files in os.walk(rootDir):
        for diri in dirs:
            if diri == 'a' and '5000ppm' in root and '333K' in root:
                #print(os.path.join(root, diri))
                plotN(os.path.join(root, diri))

directory(r"D:\Private\OneDrive - hust.edu.cn\Private\Gromacs\ILs\BmimTF2N_water_LiTF2N\BmimTFSI_sol_Litfsi")


ax[4].set_xticks(range( ddt * (tt1 // ddt), tt2, ddt))

plt.show()
fig.savefig("D:\\Private\\Desktop\\today\\" + title + ".png", dpi=600)
