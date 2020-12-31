#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 16:49:49 2020

 analysis code for processing data file

@author: yunjaeh
"""

import numpy as np
import matplotlib.pyplot as plt
# from scipy.interpolate import griddata


# %% line probes
# fPath ='./output/pcprobes/'
fPath ='./Sc_1.0/pcprobes/'
prefix = 'up.'

coord = np.loadtxt(fPath+prefix+'pxyz', skiprows=1)
nPt = len(coord)

# coord.sort(axis=0)
idx = coord[:,0]
y = coord[:,2]

u_tau = 1
nu = 1/395

y_plus=(y+1)*u_tau/nu

#%%

numSteps = 3200

U_raw=np.zeros((nPt,numSteps))
V_raw=np.zeros((nPt,numSteps))
W_raw=np.zeros((nPt,numSteps))
P_raw=np.zeros((nPt,numSteps))
CT_raw=np.zeros((nPt,numSteps))


for i in range(0,numSteps):
    print(str(i).zfill(8))
    data_temp=np.loadtxt(fPath+prefix+str(i).zfill(8)+'.pcd',skiprows=1)
    U_raw[:,i] = data_temp[:,1]
    V_raw[:,i] = data_temp[:,2]
    W_raw[:,i] = data_temp[:,3]
    P_raw[:,i] = data_temp[:,4]
    CT_raw[:,i] = data_temp[:,0]

# %%
print(U_raw.shape,idx.shape)

plt.figure(figsize=(12,4))
plt.subplot(131)
plt.plot(U_raw.mean(axis=1), y,'k.')
plt.plot(V_raw.mean(axis=1), y,'b.')
plt.plot(W_raw.mean(axis=1), y,'r.')

plt.subplot(132)
plt.plot(U_raw.std(axis=1), y,'k.')
plt.plot(V_raw.std(axis=1), y,'b.')
plt.plot(W_raw.std(axis=1), y,'r.')

plt.subplot(133)
plt.plot(CT_raw.mean(axis=1), y,'k.')
# plt.plot(CT_raw[:,3000:].mean(axis=1), y,'b.')
# plt.plot(CT_raw[:,8000:].mean(axis=1), y,'b.')
# plt.plot(CT_raw[:,10000:].mean(axis=1), y,'r.')
# pltSteps=range(0,numSteps,int(numSteps/5))
# for i in pltSteps:
    # plt.plot(CT_raw[:,i],y,'.')

plt.grid()
plt.xticks(np.arange(0,1.01,0.25))
# plt.legend(pltSteps)
    



#%%
plt.plot(y_plus,U_raw.std(axis=1),'b.')
plt.plot(y_plus,V_raw.std(axis=1),'r.')
plt.plot(y_plus,W_raw.std(axis=1),'g.')


