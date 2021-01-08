#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 17:46:49 2021

@author: yunjaeh
"""

import numpy as np
import matplotlib.pyplot as plt

means = np.loadtxt('DNS/chan395.means')
reys  = np.loadtxt('DNS/chan395.reystress')


fCase='update'
Sc='0.1'

fPath='./results_'+fCase+'/Sc_'+Sc+'/probes/'
print('Sc=', Sc, ', Output file path=', fPath)
    
coord  = np.loadtxt(fPath+'line.README')
Ux_raw = np.loadtxt(fPath+'line.u-x')
Uy_raw = np.loadtxt(fPath+'line.u-y')
Uz_raw = np.loadtxt(fPath+'line.u-z')

y  = coord[:,2]
Ux = Ux_raw[:,3:]
Uy = Uy_raw[:,3:]
Uz = Uz_raw[:,3:]


# %%
fig, ax = plt.subplots(figsize=(5,5))
ax.plot(np.mean(Ux,axis=0),y,'r')
ax.plot([means[:,2], means[:,2]], [-1+means[:,0],1-means[:,0]],'bx')
ax.set(title='Ux_mean')
ax.legend(['CharLES','DNS'])
fig.savefig('images/mean_Ux.png')

fig, [ax1,ax2,ax3] = plt.subplots(nrows=1, ncols=3,figsize=[18,4])
ax1.plot(np.std(Ux,axis=0),y,'r')
ax1.plot(np.sqrt(reys[:,2]), -1+reys[:,0],'bx')
ax1.plot(np.sqrt(reys[:,2]),  1-reys[:,0],'bx')
ax1.set(xlabel='U_rms',ylabel='y: channel height')

ax2.plot(np.std(Uy,axis=0),y,'r')
ax2.plot(np.sqrt(reys[:,3]), -1+reys[:,0],'bx')
ax2.plot(np.sqrt(reys[:,3]),  1-reys[:,0],'bx')
ax2.set(xlabel='V_rms',ylabel='y: channel height')

ax3.plot(np.std(Uz,axis=0),y,'r')
ax3.plot(np.sqrt(reys[:,4]), -1+reys[:,0],'bx')
ax3.plot(np.sqrt(reys[:,4]),  1-reys[:,0],'bx')
ax3.set(xlabel='W_rms',ylabel='y: channel height')
fig.savefig('images/rms_Uxvw.png')