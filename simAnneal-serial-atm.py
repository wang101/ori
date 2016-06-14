# -*- coding: utf-8 -*-
"""
Created on Fri May 27 12:45:39 2016

@author: wang
"""


import h5py
import numpy
import matplotlib.pyplot as plt
import pp
import os, sys, time
import pstats
import numpy as numpy
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.PDBParser import PDBParser
import pdb
import copy

def Init():
    ptask = open("task.input","r")
    para = {}
    jobs = []
    for line in ptask.readlines():
        if(line[0]=='/' or line[0]=='\n'):
            continue
        [a,b] = line.split("=")
        if a=='angle':
            jobs.append([float(x) for x in b.strip().split(',')])
        else:
            para[a]=b.strip()
    ptask.close()
    filename = para['protein_file']
    protein_name = filename.strip().split('.')[0]
    file_type = filename.strip().split('.')[1]
    if file_type == 'cif':
        mt = MMCIF2Dict(filename)
        xlist = [float(x) for x in mt['_atom_site.Cartn_x']]
        ylist = [float(x) for x in mt['_atom_site.Cartn_y']]
        zlist = [float(x) for x in mt['_atom_site.Cartn_z']]
        allarr = numpy.vstack((xlist,ylist,zlist)).T
    elif file_type == 'pdb':
        parser = PDBParser()
        structure = parser.get_structure("test", filename)
        atoms = structure.get_atoms()
        alllist = []
        xlist = []
        ylist = []
        zlist = []
        for atom in atoms:
            xlist.append(atom.get_coord()[0])
            ylist.append(atom.get_coord()[1])
            zlist.append(atom.get_coord()[2])
            alllist.append(atom.get_coord())
        allarr = numpy.array(alllist)
    if para['CENTER'] == 'ON':
        x_ave = allarr.mean(axis=0)[0]
        y_ave = allarr.mean(axis=0)[1]
        z_ave = allarr.mean(axis=0)[2]
        allarr[:,0] = allarr[:,0]-x_ave;
        allarr[:,1] = allarr[:,1]-y_ave;
        allarr[:,2] = allarr[:,2]-z_ave

    scr_size = int(para['scr_size'])
    pix_size = float(para['pix_size'])
    distance = float(para['distance'])
    wavenum = 1.0/float(para['lambda'])
    ssc = scr_size/2.0-0.5

    s = numpy.zeros((scr_size,scr_size,3))
    for i in range(scr_size):
        for j in range(scr_size):
            x = (i-ssc)*pix_size
            y = (j-ssc)*pix_size
            z = distance
            sr = numpy.sqrt(x*x+y*y+z*z)
            s[i,j,:] = numpy.array([x*wavenum/sr,y*wavenum/sr,z*wavenum/sr-wavenum])

    return s,allarr
    
def Gen_patt(s,tryarr,angles,scr_size):
    Rx = lambda t: numpy.array([[1,0,0],[0,numpy.cos(t),-numpy.sin(t)],[0,numpy.sin(t),numpy.cos(t)]])
    Ry = lambda t: numpy.array([[numpy.cos(t),0,numpy.sin(t)],[0,1,0],[-numpy.sin(t),0,numpy.cos(t)]])
    Rz = lambda t: numpy.array([[numpy.cos(t),-numpy.sin(t),0],[numpy.sin(t),numpy.cos(t),0],[0,0,1]])

    patt = numpy.zeros((scr_size,)*2,dtype=complex)
    e1,e2,e3 = angles
    r1 = Rz(e1*numpy.pi/180)
    r2 = Rx(e2*numpy.pi/180)
    r3 = Rz(e3*numpy.pi/180)
    rall = numpy.dot(numpy.dot(r1,r2),r3)
    transarr = numpy.dot(rall,tryarr.T)*1.0e-10
    patt = map(lambda x: abs(sum(numpy.exp(2*numpy.pi*1j*numpy.dot(x,transarr)))),s.reshape(scr_size*scr_size,3))
    patt = numpy.array(patt).reshape((scr_size,)*2).T
    return patt
    
def Calc_diff(patt,compare_num,patts):
    mintry = sum([numpy.linalg.norm((a[(i,)]-patt),ord='fro') for i in xrange(start,stop)])



ERR_TOL = 1.0e-3


time1 = time.time()
h = h5py.File("./lst.h5",'r')
num,scr_size = h['pattern'].shape[0:2]
patts = h['pattern'][()]
s,allarr = Init()
tryarr = copy.deepcopy(allarr)
tryarr[-1,:] = 1

iter_times = 1000
for n_iter in xrange(iter_times):
    patt = Gen_patt(s,tryarr)
    diff = Calc_diff(patt,compare_num,patts)
    tryarr