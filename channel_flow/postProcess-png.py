#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 21:42:26 2021

@author: yunjaeh
"""


import numpy as np
import matplotlib.pyplot as plt

# fCase = 'orig'
fCase = 'update'
y = dict()
CT = dict()

Sc_list = ['0.1','0.7','1.0','10.0']
for Sc in Sc_list:
    data=np.loadtxt('results_'+fCase+'/Sc_'+Sc+'/CT.mean.collapse_width.dat');
    y[Sc]  = data[:,3]
    CT[Sc] = data[:,5]
    
fig, ax = plt.subplots(figsize=(6,4))

for Sc in Sc_list:
    ax.plot(CT[Sc],y[Sc])
    
ax.set(xlabel='CT', ylabel='y: channel height')
ax.grid()
ax.legend(Sc_list)
    




