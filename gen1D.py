# -*- coding: utf-8 -*-
"""
Created on Sun May  8 10:33:33 2016

@author: wang
"""


import numpy as np
import matplotlib.pyplot as plt
import string
import os
import h5py
import copy
import pdb
import re


ptask = open("task.input","r")
job = []
para = {}
for line in ptask.readlines():
    if(line[0]=='/' or line[0]=='\n'):
        continue
    [a,b] = line.split("=")
    if a=='angle':
        job.append([float(x) for x in b.strip().split(',')])
    else:
        para[a]=b.strip()
filepath=para['save_path']
pix_size= float(para['pix_size'])
D = float(para['distance'])
lambdaxray=float(para['lambda'])
if not os.path.exists(filepath):
    os.mkdir(filepath)    
filename_list = os.listdir(filepath)
theta_num = 80
s_num = 40
s_max = 1.0/4e-10

def biInterplot(x,y,arr):
    size = arr.shape[0]
    x1 = np.floor(x)
    x2 = x1+1
    if x2 >= size:
        x2 = size-1
    y1 = np.floor(y)
    y2 = y1+1
    if y2 >= size:
        y2 = size-1
    return arr[x1,y1]*(x2-x)*(y2-y)+arr[x1,y2]*(x2-x)*(y-y1)+arr[x2,y1]*(x-x1)*(y2-y)+arr[x2,y2]*(x-x1)*(y-y1)

rere = re.compile(r'morphcto_\d+:(\w|-|\.)*:(\w|-|\.)*,(\w|-|\.)*,(\w|-|\.)*.dat')
allpatt = []
for filename in filename_list:
    if rere.match(filename.strip()):
        try:
            f = open(filepath+filename,'r')
        except:
            pass
        else:
            lst = []
            for line in f.readlines():
                arrRow = np.array([string.atof(x) for x in line.split(',')[:-1]])
                lst.append(arrRow)
            arr = np.array(lst)
            plt.imsave(filepath+filename[:-3]+".png",arr)
            
#            theta_list = [2*np.pi*x/theta_num for x in range(theta_num)]
#            s_list = [x*s_max/s_num for x in range(s_num)]
#            patt1d = []
#            mem = []
#            centerx = arr.shape[0]/2.0
#            centery = arr.shape[1]/2.0
#            for s in s_list:
#                r = np.tan(2*np.arcsin(s*lambdaxray/(4*np.pi)))*D/pix_size
#    #            pdb.set_trace()
#                sumr = 0
#                for theta in theta_list:
#                    x = r*np.cos(theta)+centerx
#                    y = r*np.sin(theta)+centery
#                    sumr = sumr + biInterplot(x,y,arr)
#                sumr=sumr*2*np.pi/theta_num;
#                patt1d.append([s,sumr])
#            patt1d = np.array(patt1d)
#            plt.plot(patt1d[:,0],patt1d[:,1])
#            plt.gca().set_xlabel("m^(-1)")
#plt.savefig(filepath+"plot1D.png")
#allpatt= np.array(allpatt)
#sumpatt = np.zeros((allpatt.shape[1],allpatt.shape[2]))
#for i in range(allpatt.shape[0]):
#    sumpatt = sumpatt + allpatt[i,:,:]
#sumpatt = sumpatt/allpatt.shape[0]
#plt.close()
#plt.plot(sumpatt[:,0],sumpatt[:,1])
#plt.savefig(filepath+"plot1Dall.png")
#plt.plot(allpatt[:,0],allpatt[:,1])