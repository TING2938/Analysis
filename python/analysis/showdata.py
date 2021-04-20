import numpy as np
from glob import glob
import matplotlib.pyplot as plt
import os

def plotN(directory):    
    os.chdir(directory)
    noLi = True
    pwd = os.getcwd().split('\\')[-3:]
    for i in 'half', 'equal', 'twice', '4v1':
        if i in pwd:
            noLi = False
            break
    
    nrow = 1 if noLi else 2
    
    fn = 'Energy-2-1*.dat', 'Energy-2-3*.dat', 'Energy-4-1*.dat', 'Energy-4-3*.dat', \
            'Energy-4-4*.dat', 'Energy-4-2*.dat', 'Energy-2-4*.dat'
    title = 'Li-Bmim', 'Li-Tf2N', 'SOL-Bmim', 'SOL-Tf2N', 'SOL-SOL', 'SOL-Li', 'Li-SOL'
    
    fn = fn[2:5] if noLi else fn
    title = title[2:5] if noLi else title
    
    fig, ax = plt.subplots(nrow, 4, figsize=(12, 3 * nrow ))
    ax = ax.flatten()
    
    ndx = 0
    for i in fn: 
        fileName = glob(i)[0]
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
        
        ax[ndx].set_xticks(range( dt * (t1 // dt), t2+1, dt))
        ax[ndx].plot(x, y, "-o", color='#FF6400', lw=1.5)
        ax[ndx].bar(x, y, color='#00A383', width = 0.200299 * np.log(2.91618 * l))
        ax[ndx].set_title(title[ndx])
        ax[ndx].set_xlabel('N (#)')
        ax[ndx].set_ylabel('percentage (%)')
        ndx += 1
        
    for i in range(ndx, 4 * nrow):
        ax[i].axis('off')
    
    Li = {'zero':'0:1', 'half':'0.5:1', 'equal':'1:1', 'twice':'2:1', '4v1':'4:1'}
    
    pwd1 = pwd.copy()
    if noLi:
        pwd = pwd[-2:]
    else:
        pwd[-2] = Li[pwd[-2]]
    fig.text(0.8, 0.2, ' / '.join(pwd))            
    fig.savefig("D:\\Private\\Desktop\\BT\\" + '_'.join(pwd1)  + '.png', dpi=600)
    plt.close()

def direct(rootDir):
    for root, dirs, files in os.walk(rootDir):
        for diri in dirs:
            if diri in ('a', 'b', 'c'):
                plotN(os.path.join(root, diri))

direct(r"D:\Private\OneDrive - hust.edu.cn\Private\Gromacs\ILs\BmimTF2N_water_LiTF2N\BmimTFSI_sol_Litfsi")