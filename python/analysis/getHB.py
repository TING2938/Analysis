# -*- coding: utf-8 -*-
import numpy as np

fn = ['HB_SOL_SOL_num.xvg', 'HB_SOL_Bmim_num.xvg', 'HB_SOL_BF4_num.xvg']

result = []

for i in fn:
    result.append(np.loadtxt(i, comments=['#', '@'])[:, 1].mean())

for i in result:
    print(f"{i:.5f}", end=' ')
