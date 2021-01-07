#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 18:20:54 2020

@author: yunjaeh
"""

#%% point cloud probe
import numpy as np

numPt = 200

X = np.zeros((numPt,))
Y = np.linspace(-1,1,numPt)
Z = np.zeros((numPt,))

coord=np.concatenate((X,Y,Z)).reshape((3,numPt))

np.savetxt('PC_line.txt',coord.transpose())

