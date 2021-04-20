# -*- coding: utf-8 -*-

def readNdx(fileName):
    """return zero-based index dict"""
    index = {}
    tmp_str = ""
    with open(fileName) as file:
        for line in file:
            if len(line) == 0 or line.isspace():
                continue
            if '[' in line and ']' in line:
                tmp_str = line[line.find('[')+1 : line.find(']')].strip()
                index[tmp_str] = []
                continue
            index[tmp_str] += [int(i)-1 for i in line.split()]
    return index

def writeNdx(fileName, data, ndxName=None, mode='a'):
    if ndxName == None:
        ndxname = fileName.split('.')[0]
    file = open(fileName, mode)
    file.write('[ ' + ndxName + ' ]\n')
    for ndx, i in enumerate(data):
        file.write("{0} ".format(i))
        if (ndx+1) % 20 == 0:
            file.write("\n")
    file.write('\n')
    file.close()
    return len(data)