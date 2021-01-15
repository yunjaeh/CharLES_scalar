#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 15:27:57 2021

@author: yunjaeh
"""

import numpy as np
import matplotlib.pyplot as plt

Pr = 0.71
mu = 0.67e-3
rho = 1
nu = mu/rho
alpha = nu/Pr
h=1

tick_list=np.arange(-0.5,0.51,0.2)
step_list = [1000,2000,3000]

y, U, T = dict(), dict(), dict()
u_rms, v_rms, w_rms = dict(), dict(), dict()

for step in step_list:
    data = np.loadtxt('out/cprobes/uvwpT.'+str(step).zfill(8)+'.cp')
    y[step] = data[:,2]
    U[step] = data[:,8]
    T[step] = data[:,20]
    u_rms[step] = data[:,5]
    v_rms[step] = data[:,9]
    w_rms[step] = data[:,13]

fig, axes = plt.subplots(nrows=2, figsize=(4,7))

for step in step_list:
    axes[0].plot(y[step],T[step])
    axes[1].plot(y[step],U[step]*h/alpha)

   
for ax in axes:
    ax.legend(step_list)
    ax.grid()
    ax.set(xticks=tick_list, xlabel='x')
axes[0].set(ylabel='T')
axes[1].set(ylabel='Uh/alpha')
    

#%% temperature and velocity profiles

step = step_list[2]
fig, axes = plt.subplots(ncols=2, figsize=(6,4))
axes[0].plot(y[step]+0.5,-T[step])
axes[1].plot(y[step]+0.5,-U[step]*h/alpha)    

for ax in axes:
    # ax.legend(step)
    ax.grid()
    ax.set(xlim=(0,0.5), xlabel='x')
axes[0].set(ylim=(0,0.5))
axes[1].set(ylim=(0,400))

# WRITE_IMAGE NAME=161065901251548 TARGET 29.392873764038086 39.36409378051758 5.315889358520508 CAMERA -127.1139144897461 180.18894958496094 342.6572265625 UP 0.16882915794849396 0.9349678158760071 -0.3119806945323944 SIZE 1856 888 WIDTH 789.6910400390625 HIDE_ZONES 0,1,3 GEOM PLANE 50 75 0 0 0 1 VAR mag(u) NO_GEOM_BLANKING GEOM ISO p 1

#%% rms values

step = step_list[2]
fig, axes = plt.subplots(figsize=(6,4))
axes.plot(y[step]+0.5,u_rms[step])
axes.plot(y[step]+0.5,v_rms[step])
axes.plot(y[step]+0.5,w_rms[step])
axes.legend(['u_rms','v_rms','w_rms'])
axes.grid()
axes.set(xlim=(0,0.5), xlabel='x', \
         ylim=(0, 0.3))





