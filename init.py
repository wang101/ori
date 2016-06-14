# -*- coding: utf-8 -*-
"""
Created on Tue May  3 16:06:13 2016

@author: wang
"""

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.PDBParser import PDBParser
import numpy as np
import matplotlib.pyplot as plt
import h5py
import os


def prob(theta,A):#theta in (-pi,pi):
   return np.exp(np.cos(theta)*A)*(np.sin(abs(theta)))

#parsing:
ptask = open("task.input","r")
para = {}
demo = 0
job = []
metropolis_loop = 500
for line in ptask.readlines():
    if(line[0]=='/' or line[0]=='\n'):
        continue
    [a,b] = line.split("=")
    if a=='angle':
        job.append([float(x) for x in b.strip().split(',')])
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
    allarr = np.vstack((xlist,ylist,zlist)).T
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
    allarr = np.array(alllist)
if para['CENTER'] == 'ON':
    x_ave = allarr.mean(axis=0)[0]
    y_ave = allarr.mean(axis=0)[1]
    z_ave = allarr.mean(axis=0)[2]
    allarr[:,0] = allarr[:,0]-x_ave;
    allarr[:,1] = allarr[:,1]-y_ave;
    allarr[:,2] = allarr[:,2]-z_ave
if para['RAND_EULER'] == 'ON':
    if para['MAG_FIELD'] == 'ON':
        mag_field = 'true'
        eulerlist = []
        eulerlist.append(np.array([0,0,0]))
        mag_x,mag_y,mag_z = para['magnetic_angle'].strip().split(',')
        mag_x = float(mag_x);mag_y = float(mag_y);mag_z = float(mag_z)
        if(mag_x*mag_y*mag_z==0):
            print "invalid magnetic direction angle, set as (0,0,1)\n"
            mag_x = 0;
            mag_y = 0;
            mag_z = 1;
        else:
            mag_all = np.sqrt(mag_x*mag_x+mag_y*mag_y+mag_z*mag_z);
            mag_x = mag_x/mag_all;
            mag_y = mag_y/mag_all;
            mag_z = mag_z/mag_all;
        temperature = float(para['temperature']);
        miu = float(para['miu']);
        tesla = float(para['tesla']);
        k =1.3806488e-23
        index = miu*tesla/(k*temperature)
        maxtheta = np.arccos((np.sqrt(1+4*index*index)-1.0)/(2*index))
        #use Rejection method sampling
        maxvalue = prob(maxtheta,index)
        valid_sample = 0
        thetalist = []
        while valid_sample < int(para['rand_euler_num']):
            proposal_theta = np.random.rand()*np.pi
            if np.random.rand() < prob(proposal_theta,index)/maxvalue:
                e1 = np.random.rand()*360#rad = radian
                e2 = proposal_theta*180/np.pi
                e3 = np.random.rand()*np.pi*2
                job.append([e1,e2,e3])
                valid_sample = valid_sample + 1
                thetalist.append(e2)
        step = 5
        xlist = np.array([x*np.pi*step/180.0 for x in range(0,int(180/step)+1)])
        ylist = np.array([prob(x,index) for x in xlist])
    else:
        mag_field = 'false'
        eulerlist = []
        eulerlist.append(np.array([0,0,0]))
        for i in range(int(para['rand_euler_num'])):
            e1 = np.random.rand()*360
            r = np.random.rand()
            e2 = np.arcsin(2*r-1)/np.pi*180
            e3 = np.random.rand()*360
            job.append([e1,e2,e3])
else:
    mag_field = 'false'
if not os.path.exists(para['save_path']):
    os.makedirs(para['save_path'])
#end parsing

#writing task file
otask = open("task","w")
otask.write(para['ground_truth']+'\n')
otask.write(para['protein_file'].split('.')[0]+'.dat\n')
otask.write(para['save_path'] + protein_name + ':\n')
otask.write(para['EULER_def'][0]+' '+para['EULER_def'][1]+' '+para['EULER_def'][2]+'\n')
for item in job:
    a,b,c = item
    otask.write(str(a)+' '+str(b)+' '+str(c)+'\n')
otask.close()
#end writing task file

#writing parameter.cpp
para4c = open("parameter.cpp","w")
para4c.write("class para\n{\
      public:\n \
           double lambda;//meter \n \
           double pix_size;//perimeter of each pixel\n \
           int demo_size; //space of 3D model\n \
           int scr_size;//screen size\n \
           double wave_number;//wavenumber = 1/lambda\n \
           double distance;//distance from crystal to screen\n \
           double epsilon;//parameter reserved for calculation of crystal\n \
           char eulerdef[3]; //def of euler angle 'zxz' as default;\n \
           int natom;//atom number\n \
           bool magnetic_field;\n\
           double magnetic_angle[3];\n\
           double tesla;\n\
           double miu;\n\
           double temperature;\n\
           para()\n \
          {\n  ")
para4c.write("\teulerdef[0]=\'{}\';\n\
            eulerdef[1]=\'{}\';\n\
            eulerdef[2]=\'{}\';\n\
            demo_size={};\n\
            scr_size={};\n\
            pix_size={};\n\
            distance={};\n\
            lambda = {};\n\
            natom = {};\n \
            magnetic_field = {};\n \
            magnetic_angle[0] = {};\n \
            magnetic_angle[1] = {};\n \
            magnetic_angle[2] = {};\n \
            tesla = {};\n \
            miu = {};\n\
            temperature = {};\n\
            ".format(\
            para['EULER_def'][0],\
            para['EULER_def'][1],\
            para['EULER_def'][2],\
            para['demo_size'],\
            para['scr_size'],\
            para['pix_size'],\
            para['distance'],\
            para['lambda'],\
            allarr.shape[0],\
            mag_field,\
            para['magnetic_angle'].strip().split(',')[0],\
            para['magnetic_angle'].strip().split(',')[1],\
            para['magnetic_angle'].strip().split(',')[2],\
            para['tesla'],\
            para['miu'],\
            para['temperature']))
para4c.write("}}")
para4c.write(";")
para4c.close()
#end writing parameter.cpp



#writing protein data(atom coordinate)ASCII
out = open(protein_name+'.dat',"w")
num = allarr.shape[0];
out.write(str(num)+'\n')
for i in range(num):
    x = allarr[i,0]
    y = allarr[i,1]
    z = allarr[i,2]
    out.write(str(x) +','+str(y) +','+str(z) +'\n')
out.close()
#end writing protein data

#writing HDF5
f = h5py.File('protein.h5','w')
jobarr = np.array(job)
f['protein_name'] = protein_name
f['euler_def'] = para['EULER_def']
f['scr_size'] = int(para['scr_size'])
f['pix_size'] = float(para['pix_size'])
f['distance'] = float(para['distance'])
f['lambda'] = float(para['lambda'])
f['CENTER'] = para['CENTER']
f['atom_num'] = num
f['job_num'] = len(job)
f.create_dataset("job_list",jobarr.shape,data=jobarr)
f.close()
#end writing HDF5


#calling c
#os.system('g++ -o main main.cpp')
#os.system('./main')
#end calling c
