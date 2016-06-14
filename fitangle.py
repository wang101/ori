# -*- coding: utf-8 -*-
"""
Created on Sun May 22 10:54:38 2016

@author: wang
"""

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
import re
import pickle as pkl

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
filename = para['protein_file']
filepath=para['save_path']
pix_size= float(para['pix_size'])
D = float(para['distance'])
lambdaxray=float(para['lambda'])
if not os.path.exists(filepath):
    os.mkdir(filepath)    
filename_list = os.listdir(filepath)

rere = re.compile(r'morphcto_\d+:(\w|-|\.)*:(\w|-|\.)*,(\w|-|\.)*,(\w|-|\.)*.dat')


arrlst = []
anglst = []
for filename in filename_list:
    if rere.match(filename.strip()):
        protein_name,filenum,angs = filename.split(":")
        ang1,ang2,ang3 = angs[:-4].split(",")
        filenum = int(filenum)
        ang1 = float(ang1)
        ang2 = float(ang2)
        ang3 = float(ang3)
        anglst.append(np.array((ang1,ang2,ang3)))
        f = open(filepath+filename,'r')
        lst = []
        for line in f.readlines():
            arrRow = np.array([string.atof(x) for x in line.split(',')[:-1]])
            lst.append(arrRow)
        arr = np.array(lst)
        arrlst.append(arr)

h = h5py.File('lst'+protein_name+'.h5','w')
sparr = arr.shape
arrall = np.array(arrlst)
angall = np.array(anglst)
protein_namelst = (protein_name,)*len(arrall)
h.create_dataset('pattern',arrall.shape,data=arrall)
h.create_dataset('angle',angall.shape,data=angall)
h.create_dataset('protein_name',(len(protein_namelst),),data=protein_namelst)
h.close()