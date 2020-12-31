#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 analysis code for processing data file

@author: yunjaeh
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


#%% line probes
fPath ='./output/probes/'
# fPath ='./Sc_0.7/probes/'
fPath='./data_cluster/Sc_0.1/output/probes/'

coord = np.loadtxt(fPath+'line.README')
y = coord[:,2]

Ux_raw = np.loadtxt(fPath+'line.u-x')
Uy_raw = np.loadtxt(fPath+'line.u-y')
Uz_raw = np.loadtxt(fPath+'line.u-z')
CT_raw = np.loadtxt(fPath+'line.CT')

Ux = Ux_raw[:,3:]
Uy = Uy_raw[:,3:]
Uz = Uz_raw[:,3:]
CT = CT_raw[:,3:]


#%% time advance of profiles
nSteps = 5000
plt.figure()
for i in range(0,len(Ux),nSteps):
    plt.subplot(121)
    plt.plot(Ux[i,:],y)
    
    plt.subplot(122)
    plt.plot(CT[i,:],y)
    # plt.plot(Uy[i,:],y)
    
    # plt.subplot(133)
    # plt.plot(Uz[i,:],y)
    # plt.plot(CT[i,:],y)

plt.subplot(121)
plt.xlabel('Ux')
plt.ylabel('y')
plt.legend(range(0,len(Ux),nSteps))


plt.subplot(122)    
plt.xlabel('CT')
plt.ylabel('y')
plt.legend(range(0,len(CT),nSteps))


# %% mean profile
init_step=15000

plt.figure(figsize=(6,4))
# plt.subplot(131)
plt.plot(y,np.mean(Ux,axis=0))
plt.plot(y,np.mean(Uy,axis=0))
plt.plot(y,np.mean(Uz,axis=0))
# plt.plot(np.mean(Ux[init_step:,],axis=0),y)
# plt.plot(np.mean(Uy[init_step:,],axis=0),y)
# plt.plot(np.mean(Uz[init_step:,],axis=0),y)
plt.xlabel('y')
plt.ylabel('U, mean')
plt.ylim(-1,21)
plt.legend(['U','V','W'])


plt.figure(figsize=(6,4))
plt.plot(y,np.std(Ux,axis=0))
plt.plot(y,np.std(Uy,axis=0))
plt.plot(y,np.std(Uz,axis=0))
# plt.plot(np.std(Ux[init_step:,],axis=0),y)
# plt.plot(np.std(Uy[init_step:,],axis=0),y)
# plt.plot(np.std(Uz[init_step:,],axis=0),y)
plt.xlabel('y')
plt.ylabel('U, rms')
plt.ylim(0,3)
plt.grid()
plt.legend(['U','V','W'])


plt.figure(figsize=(6,4))
plt.plot(y,np.mean(CT,axis=0))
plt.plot(y,np.mean(CT[init_step:,],axis=0))
plt.xlabel('y')
plt.ylabel('CT')
plt.ylim(0,1)
plt.yticks(np.arange(0,1.01,0.1))
plt.grid()



#%%
plt_step = range(0,len(CT)-1,2000)
for i in plt_step:
    plt.plot(y,np.mean(CT[i:i+2000,],axis=0))
# plt.plot(y,np.mean(CT,axis=0))
# plt.plot(y,np.mean(CT[2000:,],axis=0))
# plt.plot(y,np.mean(CT[3000:,],axis=0))

plt.legend(plt_step)
plt.xlim(-1,1)
plt.ylim(0,1)
plt.yticks(np.arange(0,1.01,0.1))
plt.grid()
plt.xlabel('y; Channel height')
plt.ylabel('CT')



#%% at specific times
fig_list = [0,10000,35000,40000]
plt.figure()
for i in fig_list:
    plt.subplot(121)
    plt.plot(Ux[i,:],y)
    
    plt.subplot(122)
    plt.plot(CT[i,:],y)
    # plt.plot(Uy[i,:],y)
    
    # plt.subplot(133)
    # plt.plot(Uz[i,:],y)
    # plt.plot(CT[i,:],y)

plt.subplot(121)
plt.xlabel('Ux')
plt.ylabel('y')
plt.legend(fig_list)


plt.subplot(122)    
plt.xlabel('CT')
plt.ylabel('y')
plt.legend(fig_list)



#%%
# fPath ='./stats/'
# fPath ='./'
# U_mean=np.loadtxt(fPath+'U.mean.raw_values.dat',comments='#')
# U_mean=np.loadtxt(fPath+'U.rms.raw_values.dat',comments='#')
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

#%%
data=np.loadtxt('/run/media/yunjaeh/Data/scratch/scalar_test/output/CT.00023000.raw_values.dat',comments='#')

X, Y = np.meshgrid(np.unique(data[:,4]), np.unique(data[:,5]))
data_grid = griddata(data[:,4:6], data[:,7],(X,Y),method='linear')

plt.figure()
# plt.fig # fig size?
# plt.subplot(121)
plt.contourf(X,Y,data_grid,100)
plt.colorbar()

plt.figure()
plt.subplot(121)
plt.plot(np.nanmean(data_grid,axis=1),np.unique(data[:,5]))
plt.grid()
# plt.xlim(0,1)

plt.subplot(122)
plt.plot(np.nanstd(data_grid,axis=1),np.unique(data[:,5]))
plt.grid()


