# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
#%%
def NR(fileName, rho, r1, r2):   
    data = np.loadtxt(fileName, comments=["#", "@"])
    r = data[:, 0]
    gr = data[:, 1]
    yy = 4 * rho * np.pi * gr * (r**2)
    nab = [np.trapz(yy[0:i], r[0:i]) for i in range(len(r))]
    
    i1 = int(r1 * len(r) / r[-1])
    i2 = int(r2 * len(r) / r[-1])
    index = np.argmin(gr[i1:i2]) + i1
    xOne, Nr = 0, 0
    for i in range(len(r)):
        if r[i] >= r1 and gr[i-1] >= 1 and gr[i] <= 1:
            xOne = (r[i-1] + r[i]) / 2
            Nr = (nab[i-1] + nab[i]) / 2
            break
    return [fileName[25:-4], r, gr, xOne, Nr, r[index], nab[index]]


def getMSD(fileName):
    result = []
    with open(fileName) as file:
        for line in file:
            if 'D[' in line:
                result.append(float(line.split()[2]))
    return result
            
            
def getLifetime(fileName):
    result = []
    with open(fileName) as file:
        for line in file:
            if "HB lifetime" in line:
                result.append(float(line.split()[-2]))
    return result
#%%
proportion = ['half', 'equal', 'twice', '4v1']
prefix = ['5000ppm/' + i + '/a/removeSOL/' for i in proportion]
fileName = ['rdf_Li_Bmim.xvg', 'rdf_Li_Tf2N.xvg']

num = 3
print('prefix: ', prefix[num])
print()

vol = np.loadtxt(prefix[num] + "volume.xvg", comments=["#", "@"])[:, 1].mean()
density = np.loadtxt(prefix[num] + "density.xvg", comments=["#", "@"])[:, 1].mean()                                         
sol = 0
li =  232

result = NR(prefix[num] + fileName[0], 500/vol, 0.8, 1.3) 
result += NR(prefix[num] + fileName[1], (500+li)/vol, 0.35, 0.8)
#%%
fig, ax = plt.subplots(2, 4, figsize=(12, 6))
ax = ax.flatten()
print("{:^10s}{:>8s}{:>8s}{:>8s}{:>8s}".format("Name", "Rone", "NRone", "Rmin", "NRmin"))
for i, j in enumerate(result):
    if i%7 == 0:
        print(f"{result[i]:^10s}{result[i+3]:8.3f}{result[i+4]:8.3f}{result[i+5]:8.3f}{result[i+6]:8.3f}")
        ax[i//7].hlines(1, 0, 2, 'r', ':')
        ax[i//7].plot(result[i+1], result[i+2], label=j)
        ax[i//7].legend(frameon=False)
ax[7].axis("off")
fig.show()

print(f"\n{density:.3f}", end=' ')
    

print(f"{vol:.3f}", end=' ')


for i in range(len(result)): 
    if i%7 in (3, 4, 5, 6) : print(f"{result[i]:.3f}", end=' ')
    