#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 11:15:34 2021

@author: yunjaeh
"""

import numpy as np
import matplotlib.pyplot as plt

Pr      = 0.71
rho     = 1
h       = 1
dT      = 1
g       = 10
beta    = 0.0034

mu      = np.array([2.11e-4, 1.0986e-4, 6.948e-5])
nu      = mu/rho
alpha   = nu/Pr    
Rayleigh= g*beta*dT*h**3.0 / (nu*alpha)
Ra_case = ['5.4e5','2.0e6','5.0e6']

#%%
tick_list=np.arange(-0.5,0.51,0.1)

step=27000

fPath = '/home/yunjaeh/Codes/CharLES_scalar/vertical_channel/'
y, U, T = dict(), dict(), dict()
u_rms, v_rms, w_rms = dict(), dict(), dict()

for Ra in Ra_case:
    print(Ra)
    data = np.loadtxt(fPath+'Ra_'+Ra+
                      '/output/cprobes/uvwpT.'+str(step).zfill(8)+'.cp')
    y[(Ra,step)] = data[:,2]
    U[(Ra,step)] = data[:,8]
    T[(Ra,step)] = data[:,20]
    u_rms[(Ra,step)] = data[:,5]
    v_rms[(Ra,step)] = data[:,9]
    w_rms[(Ra,step)] = data[:,13]

#%% plot 
ticks_channel=[-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5]

# temperature profile
fig, ax = plt.subplots(nrows=1, figsize=(5,4))
for Ra in Ra_case:
    ax.plot(y[(Ra,step)],T[(Ra,step)])

ax.legend(Ra_case)
ax.grid()
ax.set(xlim=(-0.5,0.5), xlabel='x', xticks=ticks_channel,
       ylim=(-0.5,0.5), ylabel='T', yticks=ticks_channel)

# velocity profile
fig, ax = plt.subplots(nrows=1, figsize=(5,4))    
for i, Ra in enumerate(Ra_case):
    print(i, Ra)
    ax.plot(y[(Ra,step)],U[(Ra,step)]*h/alpha[i],[i])

    
ax.legend(Ra_case)
ax.grid()
ax.set(xlabel='x', ylabel='U h / alpha', 
       xticks=ticks_channel, yticks=np.arange(-1000,1001,100), 
       xlim=(-0.5,0.5), ylim=(-1000,1000))


#%% half profiles 

fig, ax = plt.subplots(nrows=1, figsize=(5,4))
cid='rbk'
for i, Ra in enumerate(Ra_case):
    ax.plot(y[(Ra,step)],-T[(Ra,step)],cid[i])
for i, Ra in enumerate(Ra_case):    
    ax.plot(-y[(Ra,step)]+0.5,T[(Ra,step)],cid[i])
ax.set(xlabel='x', ylabel='T',
       xlim=(0,0.5), ylim=(0,0.5))
ax.grid()
ax.legend(Ra_case)

cid='rbk'
ig, ax = plt.subplots(nrows=1, figsize=(5,4))
for i, Ra in enumerate(Ra_case):
    print(i, Ra)
    ax.plot(y[(Ra,step)]+0.5,-U[(Ra,step)], cid[i])
    # ax.plot(y[(Ra,step)]+0.5,-U[(Ra,step)]*h/alpha[i], cid[i])
for i, Ra in enumerate(Ra_case):
    ax.plot(-y[(Ra,step)]+0.5,U[(Ra,step)], cid[i])
    # ax.plot(-y[(Ra,step)]+0.5,U[(Ra,step)]*h/alpha[i], cid[i])

ax.set(xlabel='x', ylabel='U',
        xlim=(0,0.5) , ylim=(0,0.1))
ax.set(xlabel='x', ylabel='U h /alpha ',
       # xlim=(0,0.5) , ylim=(0,1000), yticks=range(0,1001,100))
ax.grid()
ax.legend(Ra_case)

    

# %% temperature and velocity half profiles

cid='rbk'
fig, ax = plt.subplots(figsize=(6,4))
for i, Ra in enumerate(Ra_case):
    ax.plot(y[(Ra,step)]+0.5,(u_rms[(Ra,step)]*h/alpha[i])**2/float(Ra)**(8.0/9.0),
            cid[i]+'-')
    ax.plot(y[(Ra,step)]+0.5,(v_rms[(Ra,step)]*h/alpha[i])**2/float(Ra)**(8.0/9.0),
            cid[i]+'--')
    ax.plot(y[(Ra,step)]+0.5,(w_rms[(Ra,step)]*h/alpha[i])**2/float(Ra)**(8.0/9.0),
            cid[i]+'.-')
ax.legend(['uu','vv','ww'])
ax.grid()
ax.set(xlim=(0,0.5), xlabel='x', \
       ylim=(0,0.6), ylabel='Reynolds stress')

    

#%%

step = step_list[-2]
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(y[step],-T[step])
ax.grid()

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




