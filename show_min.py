# -*- coding: utf-8 -*-
"""
Created on Fri May 27 12:45:39 2016

@author: wang
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

num = 10000
h = h5py.File("lst.h5",'r')
#std = h['3965:275.031,14.3841,0.dat']

std = 9828
name = h['prop'+str(std)][()]
theta = float(name[:-4].split(':')[1].split(',')[0])*np.pi/180.0
fai = float(name[:-4].split(':')[1].split(',')[1])*np.pi/180.0
sx = np.cos(fai)*np.sin(theta)
sy = np.cos(fai)*np.cos(theta)
sz = np.sin(fai)

ldist = []
ldiff = []
lx =[]
ly =[]
lz = []
for i in range(10000):
    name = h['prop'+str(i)][()]
    theta = float(name[:-4].split(':')[1].split(',')[0])*np.pi/180.0
    fai = float(name[:-4].split(':')[1].split(',')[1])*np.pi/180.0
    x = np.cos(fai)*np.sin(theta)
    y = np.cos(fai)*np.cos(theta)
    z = np.sin(fai)
    dist = np.sqrt((sx-x)**2+(sy-y)**2+(sz-z)**2)
    dis = 2*np.arcsin(dist/2)/np.pi*180
    if 1:
        ldist.append(dis)
        lx.append(x)
        ly.append(y)
        lz.append(z)
        ldiff.append(h['diff'+str(i)][()])
#    
plt.scatter(ldist,ldiff)
a = plt.gca()
a.set_xlim([0,180])
        
#fig = plt.figure()
#ax = fig.add_subplot(111,projection='3d')
#ax.scatter(lx,ly,lz,c=ldiff)


#plt.hist2d(ldist,ldiff,bins=40)