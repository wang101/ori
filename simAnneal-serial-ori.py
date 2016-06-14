# -*- coding: utf-8 -*-
"""
Created on Fri May 27 12:45:39 2016

@author: wang
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import pp
import os, sys, time
import pstats

def find_nbr(ltheta,lfai,start,stop,newt,newf):
    #use L1 distance to replace L2 distance, make 2x faster
#    following 3 lines of code are slow.
#    L1distance = lambda theta,fai:abs((theta-newt))+abs((fai-newf))
#    val = map(L1distance,ltheta,lfai)
#    mini = np.argmin(val)
    mini = -1
    minval = 100000
    for i in xrange(start,stop):
        theta = ltheta[i]
        fai = lfai[i]
        mintry = abs((theta-newt))+abs((fai-newf))
        if mintry < minval:        
            minval = mintry
            mini = i
    return mini
    

num = 10000
h = h5py.File("lst.h5",'r')
#std = h['3965:275.031,14.3841,0.dat']
time1 = time.time()

std = 9828
name = h['prop'+str(std)][()]
theta = float(name[:-4].split(':')[1].split(',')[0])*np.pi/180.0
fai = float(name[:-4].split(':')[1].split(',')[1])*np.pi/180.0
sx = np.cos(fai)*np.sin(theta)
sy = np.cos(fai)*np.cos(theta)
sz = np.sin(fai)
sarr = h['arr'+str(std)][:]

ltheta = []
lfai = []
lname = []
ldiff = []
for i in xrange(10000):
    name = h['prop'+str(i)][()]
    lname.append(name)
    diff = h['diff'+str(i)][()]
    ldiff.append(diff)
    theta = float(name[:-4].split(':')[1].split(',')[0])*np.pi/180.0
    ltheta.append(theta)
    fai = float(name[:-4].split(':')[1].split(',')[1])*np.pi/180.0
    lfai.append(fai)
tname = tuple(lname)
ttheta = tuple(ltheta)
tfai = tuple(lfai)

lx =[]
ly =[]
lz = []
nowt = 0#now theta
nowf = 0#now fai
nowdiff = 1000000
slogdis = []
slogdiff =[]
slogrand = []
slogtemp = []
slogstep = []
locmin = 0

temp = 20
step = 1

print str(time.time()-time1)

for n_iter in xrange(5000):
    newt = nowt + (np.random.rand()-0.5)*step
    newf = nowf + (np.random.rand()-0.5)*step
    if newt<0:
        newt += np.pi*2
    elif newt>np.pi*2:
        newt -= np.pi*2
    if newf<-np.pi:
        newf += np.pi*2
    elif newf>np.pi:
        newf -= np.pi*2
    minname = ""

    start = 0
    stop = 10000
    mini = find_nbr(ltheta,lfai,start,stop,newt,newf)
    
    name = lname[mini]
    theta = ltheta[mini]
    fai = lfai[mini]
    diff = ldiff[mini]
    
    x = np.cos(fai)*np.sin(theta)
    y = np.cos(fai)*np.cos(theta)
    z = np.sin(fai)
    dist = np.sqrt((sx-x)**2+(sy-y)**2+(sz-z)**2)
    dis = 2*np.arcsin(dist/2)
    
    rand_num = np.random.rand()
    if diff<nowdiff:
        temp = temp*(diff/nowdiff)**(1/4.0)
        slogtemp.append(temp)
        step = step*(diff/nowdiff)**(1/8.0)
        slogstep.append(step)
    if diff<nowdiff or rand_num < np.exp(-(diff-nowdiff)/temp):
        nowdiff = diff
        nowt = newt
        nowf = newf 
        slogdis.append(dis)
        slogdiff.append(diff)
        slogrand.append(rand_num+1*(nowdiff>diff))
        if len(slogdiff)>10 and diff > 100 and np.std(slogdiff[-7:])<7 :
            temp = 10
            locmin +=1
            step = 1
    if temp == 0:
        break
    if np.mod(n_iter,500)==0:
        print str(n_iter)+'\t'+str(diff)+'\t' +str(temp)

plt.subplot(121)
plt.plot(range(len(slogdis)),slogdis)
plt.subplot(122)
plt.plot(range(len(slogdiff)),slogdiff)
#fig = plt.figure()
#ax = fig.add_subplot(111,projection='3d')
#ax.scatter(lx,ly,lz,c=ldiff)

dur = time.time()-time1
print 'duration\tn_iter\tpertime'
print str(dur)  + '\t' +str(n_iter)+'\t'+str(dur/n_iter)

#plt.hist2d(ldist,ldiff,bins=40)