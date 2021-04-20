# -*- coding: utf-8 -*-
import numpy as np
import time

class ReadGro:
    def __init__(self, fileName):
        self.fileName = fileName
        file = open(fileName)
        self.title = file.readline().strip()
        self.natoms = int(file.readline())
        self.resnum = np.zeros(self.natoms, dtype=int)
        self.resname = np.zeros(self.natoms, dtype='U5')
        self.atomname = np.zeros(self.natoms, dtype='U5')
        self.x = np.zeros((self.natoms, 3))
        for i in range(self.natoms):
            line = file.readline()
            self.resnum[i] = int(line[0:5])
            self.resname[i] = line[5:10].strip()
            self.atomname[i] = line[10:15].strip()
            self.x[i, 0] = float(line[20:28])
            self.x[i, 1] = float(line[28:36])
            self.x[i, 2] = float(line[36:44])
        self.box = np.loadtxt(file)
        file.close()
        self.allLine = open(fileName).readlines()

    def save(self, fileName):
        file = open(fileName, 'w', newline='\n')
        file.write("{:s}\n".format(self.title))
        file.write("{:>5d}\n".format(self.natoms))
        for i in range(self.natoms):
            file.write("{0:5d}{1:<5s}{2:>5s}{3:5d}{4[0]:8.3f}{4[1]:8.3f}{4[2]:8.3f}\n"
                       .format(self.resnum[i], self.resname[i], self.atomname[i],
                               i+1, self.x[i]))
        file.write("{0[0]:10.5f}{0[1]:10.5f}{0[2]:10.5f}\n".format(self.box))
        file.close()


class WriteGro:
    def __init__(self, fnm):
        self.__file = open(fnm, 'w', newline='\n')
        self.__file.write('Created by python(itp.gro). {:s}\n'.format(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime())))
        self.__filePos = self.__file.tell()
        self.__file.write('{:5s}\n'.format(' '))
        self.natoms = 0

    def write(self, data, resname, atomname, resnum="auto"):
        if resnum == "auto":
            resnum = list(range(self.natoms+1, self.natoms+len(data)+1))
        elif isinstance(resnum, int):
            resnum = np.full(len(data), resnum)
        if isinstance(resname, str):
            resname = np.full(len(data), resname)
        if isinstance(atomname, str):
            atomname = np.full(len(data), atomname)

        for i in range(len(data)):
            self.natoms += 1
            self.__file.write("{:5d}{:<5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n"
                            .format(resnum[i], resname[i], atomname[i], self.natoms,
                                    data[i, 0], data[i, 1], data[i, 2]))

    def finish(self, box):
        self.box = box
        self.__file.write("{0[0]:10.5f}{0[1]:10.5f}{0[2]:10.5f}\n".format(box))
        self.__file.seek(self.__filePos)
        self.__file.write('{:5d}\n'.format(self.natoms))
        self.__file.close()
