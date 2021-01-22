#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 15:27:57 2021

@author: yunjaeh
"""

import numpy as np
import matplotlib.pyplot as plt

Pr      = 0.71
# mu      = 6.7e-4
# mu      = 2.11e-4     # Ra = 5.4e5
# mu      = 1.0986e-4   # Ra = 2.0e6
mu      = 6.948e-5  # Ra = 5.0e6
rho     = 1
nu      = mu/rho
alpha   = nu/Pr
h       = 1
dT      = 1
g       = 10
beta    = 0.0034
Ra      = g*beta*dT*h**3.0 / (nu*alpha)

print(Ra)
tick_list=np.arange(-0.5,0.51,0.1)
# step_list = [2000,5000,6000,10000,11000,12000,13000]
step_list = [0,21000,23000,25000,27000]
# step_list = [0,10000,11000,15000]
# step_list = range(0,9001,2000)

y, U, T = dict(), dict(), dict()
u_rms, v_rms, w_rms = dict(), dict(), dict()

for step in step_list:
    data = np.loadtxt('cprobes/uvwpT.'+str(step).zfill(8)+'.cp')
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
    

 #%% temperature and velocity half profiles
ticks_channel=[-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5]

step = step_list[-2]
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(y[step],-T[step])
ax.grid()
ax.set(xlim=(-0.5,0.5), xlabel='x', xticks=ticks_channel,
       ylim=(-0.5,0.5), ylabel='T', yticks=ticks_channel)
# fig.savefig('results/T_profile.png')

fig, ax = plt.subplots(figsize=(6,4))
ax.plot(y[step],-U[step]*h/alpha)    
ax.grid()
ax.set(xlim=(-0.5,0.5), xlabel='x',xticks=ticks_channel,
       ylim=(-500,500), ylabel='U')
# fig.savefig('results/U_profile.png')

fig, axes = plt.subplots(figsize=(6,4))
axes.plot(y[step]+0.5,u_rms[step]*h/alpha)
axes.plot(y[step]+0.5,v_rms[step]*h/alpha)
axes.plot(y[step]+0.5,w_rms[step]*h/alpha)
axes.legend(['u_rms','v_rms','w_rms'])
axes.grid()
axes.set(xlim=(0,0.5), xlabel='x', xticks=ticks_channel[6:-1],
         ylim=(0,300), ylabel='U rms')
# fig.savefig('results/U_rms.png')

#%% temperature and velocity half profiles

# step = step_list[-1]
# fig, ax = plt.subplots(figsize=(6,4))
# ax.plot(y[step]+0.5,-T[step])
# ax.grid()
# ax.set(xlim=(0,0.5), xlabel='x', 
#        ylim=(0,0.5), ylabel='T')

# fig, ax = plt.subplots(figsize=(6,4))
# ax.plot(y[step]+0.5,-U[step]*h/alpha)    
# ax.grid()
# ax.set(xlim=(0,0.5), xlabel='x',
#        ylim=(0,400), ylabel='U')


fig, axes = plt.subplots(figsize=(6,4))
axes.plot(y[step]+0.5,(u_rms[step]*h/alpha)**2.0/Ra**(8.0/9.0))
axes.plot(y[step]+0.5,(v_rms[step]*h/alpha)**2.0/Ra**(8.0/9.0))
axes.plot(y[step]+0.5,(w_rms[step]*h/alpha)**2.0/Ra**(8.0/9.0))
axes.legend(['u_rms','v_rms','w_rms'])
axes.grid()
axes.set(xlim=(0,0.5), xlabel='x', \
         ylim=(0,0.6), ylabel='U rms')


# # Reynolds stresses
# scale_factor = Ra **(8.0/9.0)
# fig, axes = plt.subplots(figsize=(6,4))
# axes.plot(y[step]+0.5,(u_rms[step]*h/alpha)**2.0/scale_factor)
# axes.plot(y[step]+0.5,(v_rms[step]*h/alpha)**2.0/scale_factor)
# axes.plot(y[step]+0.5,(w_rms[step]*h/alpha)**2.0/scale_factor)
# axes.legend(['u''u''','v''v''','w''w'''])
# axes.grid()
# axes.set(xlim=(0,0.5), xlabel='x', \
#          ylim=(0,0.45), ylabel='Reynolds stress')




