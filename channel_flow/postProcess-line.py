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
Sc_list = ['0.1', '0.7', '1.0', '10.0']

y, Ux, Uy, Uz, CT =dict(), dict(), dict(), dict(), dict()

for Sc in Sc_list:
    fPath='./output/Sc_'+Sc+'/probes/'
    print('Sc=', Sc, ', Output file path=', fPath)
    
    coord  = np.loadtxt(fPath+'line.README')
    Ux_raw = np.loadtxt(fPath+'line.u-x')
    Uy_raw = np.loadtxt(fPath+'line.u-y')
    Uz_raw = np.loadtxt(fPath+'line.u-z')
    CT_raw = np.loadtxt(fPath+'line.CT')

    y[Sc]  = coord[:,2]
    Ux[Sc] = Ux_raw[:,3:]
    Uy[Sc] = Uy_raw[:,3:]
    Uz[Sc] = Uz_raw[:,3:]
    CT[Sc] = CT_raw[:,3:]

# %% mean profile
init_step=10000

# velocity
Sc='0.7'
plt.figure(1,figsize=(5,5))
plt.plot(np.mean(Ux[Sc][init_step:,:],axis=0),y[Sc])
plt.plot(np.mean(Uy[Sc][init_step:,:],axis=0),y[Sc])
plt.plot(np.mean(Uz[Sc][init_step:,:],axis=0),y[Sc])
plt.title('Mean velocity profile: Sc='+Sc)
plt.legend(['Ux','Uy','Uz'])
plt.xlabel('U: velocity')
plt.ylabel('y: channel height')
plt.grid()
plt.savefig('./images/velocity_profile_mean.png')

plt.figure(2,figsize=(5,5))
plt.plot(np.std(Ux[Sc][init_step:,:],axis=0),y[Sc])
plt.plot(np.std(Uy[Sc][init_step:,:],axis=0),y[Sc])
plt.plot(np.std(Uz[Sc][init_step:,:],axis=0),y[Sc])
plt.title('RMS velocity profile: Sc='+Sc)
plt.legend(['Ux','Uy','Uz'])
plt.xlabel('U: velocity')
plt.ylabel('y: channel height')
plt.grid()
plt.savefig('images/velocity_profile_rms.png')


#%% Scalar 

plt.figure(3,figsize=(5,5))
for Sc in Sc_list:
    plt.plot(np.mean(CT[Sc][init_step:,:],axis=0),y[Sc])

plt.title('Mean scalar profile')
plt.xlabel('CT: scalar')
plt.ylabel('y: channel height')
plt.xticks(np.arange(0,1.01,0.1))
plt.grid()
plt.legend(Sc_list)
plt.savefig('images/scalar_profile_mean.png')

 
#%% time advance of profiles

Sc='0.7'
nStep=5000
plt_step = range(0,len(CT[Sc])-1,nStep)

plt.figure(4,figsize=(5,5))
for i in plt_step:
    plt.plot(np.mean(CT[Sc][i:i+nStep,],axis=0),y[Sc])
plt.title('Time advance of scalar profile, Sc='+Sc)
plt.xlim(0,1)
plt.ylim(-1,1)
plt.xlabel('CT')
plt.ylabel('y; Channel height')
plt.xticks(np.arange(0,1.01,0.1))
plt.legend(plt_step)
plt.grid()
plt.savefig('images/scalar_time_advance.png')


