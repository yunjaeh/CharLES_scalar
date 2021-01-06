#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 17:38:55 2021

@author: yunjaeh
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

ticks=[0, 0.5, 1]
T_max = 100

for i in ['4','5','6']:
# for i in ['5']:
    print('Process data: Ra=10^'+i)
    dataT=np.loadtxt('results/T.Ra_10_'+i+'.dat')
    dataU=np.loadtxt('results/U.Ra_10_'+i+'.dat')
    dataV=np.loadtxt('results/V.Ra_10_'+i+'.dat')
    
    # 1. non-uniform grid
    # X, Y = np.meshgrid(np.unique(data[:,4]), np.unique(data[:,5]))
    # 2. uniform grid
    X, Y = np.meshgrid(np.arange(0,1,0.01),np.arange(0,1,0.01))

    # griddata for contour plots & streamlines
    T = griddata(dataT[:,4:6], dataT[:,7],(X,Y),method='linear')
    U = griddata(dataU[:,4:6], dataU[:,7],(X,Y),method='linear')
    V = griddata(dataV[:,4:6], dataV[:,7],(X,Y),method='linear')

    fig, [ax1,ax2] = plt.subplots(nrows=1,ncols=2, figsize=(12,5))

    # streamline
    ax1.streamplot(X,Y,U,V,density=2,minlength=0.5)
    ax1.set(title='Streamline',xlim=(0,1),ylim=(0,1),xticks=ticks,yticks=ticks)

    # Temperature contour
    cf=ax2.contourf(X,Y,T/T_max,100,cmap=plt.cm.jet)
    fig.colorbar(cf,ax=ax2,ticks=[0,0.5,1])
    ax2.contour(cf,colors='k',levels=np.arange(0,1,0.1))
    ax2.set(title='Temperature, Ra=10^'+i,xticks=ticks,yticks=ticks)
    
    fig.savefig('images/results_Ra_10_'+i+'.png')









