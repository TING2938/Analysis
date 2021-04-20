#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 18:12:16 2019
@author: TING 
"""
import numpy as np
import matplotlib.pyplot as plt


#%%
# parameter:
L = 1
R = 1.139
C = 9.922
Mass = 0.0123

u0 = 4.0
P = 4.0

T = 1
n = 1000
m = 1000

h = L / n
tao = T / m
a = 1 / (R * C)
#%%

u = np.zeros((n, m))

# 1.
u[:, 0] = u0
u[:, 1] = u0
u[0, 1] = np.sqrt(u[1, 1]**2 - 2 * h * R * P)

for j in range(2, m-10):
    for i in range(1, n-1):
        u[i, j] = 2 * tao * a / h / h * (u[i+1, j-1] - 2 * u[i, j-1] + u[i-1, j-1]) + u[i, j-2]
    u[0, j] = np.sqrt(u[1, j]**2 - 2 * h * R * P)
    u[n-1, j] = u[n-2, j]


