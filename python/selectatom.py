# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 20:19:10 2020

@author: TING
"""

resname = 'H'

newfile = open('new.gro', 'w')
file = open('nvt.gro')

newfile.write(file.readline())
natoms = int(file.readline())
newfile.write("{:5d}\n".format(natoms))

for i in range(natoms):
    line = file.readline()
    if line[5:10].strip() != resname:
        newfile.write(line)
    
import matplotlib.pyplot as plt


            
