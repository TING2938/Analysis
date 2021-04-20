# -*- coding: utf-8 -*-
import ctypes
import time
import numpy as np
from numpy.ctypeslib import ndpointer
dll = ctypes.windll.LoadLibrary("libfunc/x64/Debug/libfunc.dll")
#%% primeQ
dll.primeQ.argtypes = [ctypes.c_ulonglong]
dll.primeQ.restype = ctypes.c_ulonglong
primeQ = dll.primeQ
print(primeQ(105))
#%%
dll.speed.argtypes = [ctypes.c_ulonglong]
dll.speed.restype = ctypes.c_ulonglong
#%%
beg = time.time()
print(dll.speed(5179357973))
print(time.time() - beg)

#%%
from speed import speed

import time
beg = time.time()

print(speed(5179357973))

print(time.time() - beg)
#%%
SUM = dll.SUM
SUM.argtypes = [ctypes.c_ulonglong, ctypes.c_ulonglong]
SUM.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_double))

