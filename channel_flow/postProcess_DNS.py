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
Sc='0.7'

fPath='./output/'+fCase+'/Sc_'+Sc+'/output/images_stats/'
print('Sc=', Sc, ', Output file path=', fPath)
    
U_mean = np.loadtxt(fPath+'U.mean.collapse_width.dat')
U_rms = np.loadtxt(fPath+'U.rms.collapse_width.dat')
# V_rms = np.loadtxt(fPath+'V.rms.collapse_width.dat')
# W_rms = np.loadtxt(fPath+'W.rms.collapse_width.dat')

fig, [ax1,ax2] = plt.subplots(ncols=2, figsize=(12,5))
ax1.plot(U_mean[:,3],U_mean[:,5],'r')
ax1.plot([-1+means[:,0],1-means[:,0]],[means[:,2], means[:,2]],'bx')
ax1.set(xlabel='y: channel height', ylabel='U_mean')
ax1.legend(['CharLES','DNS'])


ax2.plot(U_rms[:,3],U_rms[:,5],'r')
ax2.plot([-1+reys[:,0],1-reys[:,0]],np.sqrt([reys[:,2], reys[:,2]]),'bx')
ax2.set(xlabel='y: channel height', ylabel='U_rms')
fig.savefig('images/vel_mean_rms.png')


''' for rms v, w comp
fig, [ax1,ax2,ax3] = plt.subplots(nrows=1, ncols=3,figsize=[18,4])
ax1.plot(U_rms[:,5],U_rms[:,3],'r')
ax1.plot(U_rms[:,6],U_rms[:,3],'r')
ax1.plot(U_rms[:,7],U_rms[:,3],'r')
ax1.plot(np.sqrt(reys[:,2]), -1+reys[:,0],'bx')
ax1.plot(np.sqrt(reys[:,2]),  1-reys[:,0],'bx')
ax1.set(xlabel='U_rms',ylabel='y: channel height')

ax2.plot(V_rms[:,5],V_rms[:,3],'r')
ax2.plot(V_rms[:,6],V_rms[:,3],'r')
ax2.plot(V_rms[:,7],V_rms[:,3],'r')
ax2.plot(np.sqrt(reys[:,3]), -1+reys[:,0],'bx')
ax2.plot(np.sqrt(reys[:,3]),  1-reys[:,0],'bx')
ax2.set(xlabel='V_rms',ylabel='y: channel height')

ax3.plot(W_rms[:,5],W_rms[:,3],'r')
ax3.plot(W_rms[:,6],W_rms[:,3],'r')
ax3.plot(W_rms[:,7],W_rms[:,3],'r')
ax3.plot(np.sqrt(reys[:,4]), -1+reys[:,0],'bx')
ax3.plot(np.sqrt(reys[:,4]),  1-reys[:,0],'bx')
ax3.set(xlabel='W_rms',ylabel='y: channel height')
# fig.savefig('images/rms_Uxvw.png')
'''