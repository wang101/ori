# -*- coding: utf-8 -*-
"""
Created on Fri May 27 12:45:39 2016

@author: wang
June 2nd:
We can improve this by:
1) applying dynamic number of cores
2) combine simulation annealing and sto grad dec

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

#global slogdis,slogdiff, slogrand, slogtemp, slogstep, sloglocmin
#
#slogdis = []
#slogdiff =[]
#slogrand = []
#slogtemp = []
#slogstep = []
#sloglocmin = 0

def init():
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


def find_nbr(start,stop,s,allarr,newt,newf,new_crd):
    #use L1 distance to replace L2 distance, make 2x faster
#    following 3 lines of code are slow.
#    L1distance = lambda theta,fai:abs((theta-newt))+abs((fai-newf))
#    val = map(L1distance,ltheta,lfai)
#    mini = numpy.argmin(val)
    mini = -1
    minval = 100000
    scr_size = 64

    Rx = lambda t: numpy.array([[1,0,0],[0,numpy.cos(t),-numpy.sin(t)],[0,numpy.sin(t),numpy.cos(t)]])
    Ry = lambda t: numpy.array([[numpy.cos(t),0,numpy.sin(t)],[0,1,0],[-numpy.sin(t),0,numpy.cos(t)]])
    Rz = lambda t: numpy.array([[numpy.cos(t),-numpy.sin(t),0],[numpy.sin(t),numpy.cos(t),0],[0,0,1]])

    patt = numpy.zeros((scr_size,)*2,dtype=complex)
    e1,e2,e3 = new_crd
    r1 = Rz(e1*numpy.pi/180)
    r2 = Rx(e2*numpy.pi/180)
    r3 = Rz(e3*numpy.pi/180)
    rall = numpy.dot(numpy.dot(r1,r2),r3)
    transarr = numpy.dot(rall,allarr.T)*1.0e-10
    patt = map(lambda x: abs(sum(numpy.exp(2*numpy.pi*1j*numpy.dot(x,transarr)))),s.reshape(scr_size*scr_size,3))
    patt = numpy.array(patt).reshape((scr_size,)*2).T

    h = h5py.File("./lst.h5",'r')
    a = h['pattern'][()]
    for i in xrange(start,stop):
        patt_exp = a[(i,)]
#        patt_exp = h['arr'+str(i)]
        mintry = numpy.linalg.norm((patt_exp-patt),ord='fro')
        if mintry < minval:
            minval = mintry
            mini = i
    return mini,minval

def one_core(itertime,ltheta,lfai,nowt,nowf,nowdiff,temp,step,s,allarr,sx,sy,sz,now_crd):
    slogdis = []
    slogdiff =[]
    slogrand = []
    slogtemp = []
    slogstep = []
    sloglocmin = 0
    now1,now2,now3 = now_crd
    for n_iter in range(itertime):
        newt = nowt + (numpy.random.rand()-0.5)*step
        newf = nowf + (numpy.random.rand()-0.5)*step
        new1 = now1 + (numpy.random.rand()-0.5)*step
        new2 = now2 + (numpy.random.rand()-0.5)*step
        new3 = now3 + (numpy.random.rand()-0.5)*step
        new_crd = (new1,new2,new3)
        if newt<0:
            newt += numpy.pi*2
        elif newt>numpy.pi*2:
            newt -= numpy.pi*2
        if newf<-numpy.pi:
            newf += numpy.pi*2
        elif newf>numpy.pi:
            newf -= numpy.pi*2
        start = 0
        stop = 10000
        mini,diff = find_nbr(start,stop,s,allarr,newt,newf,new_crd)

        fai = lfai[mini]
        theta = ltheta[mini]

        x = numpy.cos(fai)*numpy.sin(theta)
        y = numpy.cos(fai)*numpy.cos(theta)
        z = numpy.sin(fai)
        dist = numpy.sqrt((sx-x)**2+(sy-y)**2+(sz-z)**2)
        dis = 2*numpy.arcsin(dist/2)

        rand_num = numpy.random.rand()
        if diff<nowdiff:
            temp = temp*(diff/nowdiff)**(1/4.0)
            slogtemp.append(temp)
            step = step*(diff/nowdiff)**(1/8.0)
            slogstep.append(step)
        if diff<nowdiff or rand_num < numpy.exp(-(diff-nowdiff)/temp):
            nowdiff = diff
            nowt = newt
            nowf = newf
            slogdis.append(dis)
            slogdiff.append(diff)
            slogrand.append(rand_num+1*(nowdiff>diff))
            if len(slogdiff)>10 and diff > 100 and numpy.std(slogdiff[-7:])<7 :
                temp = 10
#                sloglocmin +=1
                step = 1
        if temp == 0:
            break
#        if numpy.mod(n_iter,500)==0:
#            print str(n_iter)+'\t'+str(diff)+'\t' +str(temp)
    return (nowdiff,nowt,nowf,temp,step,slogdiff,slogdis,n_iter)


num = 10000
h = h5py.File("./lst.h5",'r')
#std = h['3965:275.031,14.3841,0.dat']
time1 = time.time()

std = 9828
name = h['prop'+str(std)][()]
theta = float(name[:-4].split(':')[1].split(',')[0])*numpy.pi/180.0
fai = float(name[:-4].split(':')[1].split(',')[1])*numpy.pi/180.0
print 'std: t:' + str(theta) + '\tf:' + str(fai)
sx = numpy.cos(fai)*numpy.sin(theta)
sy = numpy.cos(fai)*numpy.cos(theta)
sz = numpy.sin(fai)
#sarr = h['arr'+str(std)][:]

ltheta = []
lfai = []
for i in xrange(10000):
    theta = h['angle'][(i,0)]*numpy.pi/180.0
    ltheta.append(theta)
    fai = h['angle'][(i,1)]*numpy.pi/180.0
    lfai.append(fai)

lx =[]
ly =[]
lz = []
nowt = 0#now theta
nowf = 0#now fai
nowdiff = 1000000
temp = 20
step = 1
s,allarr = init()

print 'init time:' + str(time.time()-time1)



ppservers = ()
if len(sys.argv) > 1:
    ncpus = int(sys.argv[1])
else:
    ncpus = 3
job_server = pp.Server(ncpus, ppservers=ppservers)
jobs = []

all_n_iter = 5000
check_iter = 500
all_real_iter = 0
glogdiff = []
glogdis = []
glog = []
jobs = [0]*ncpus

now_crd = (5.0,7.0,7.0)
for n_iter in range(50):
    nowdiff,nowt,nowf,temp,step,ldiff,ldis,n_re_iter = one_core(1,ltheta,lfai,nowt,nowf,nowdiff,temp,step,s,allarr,sx,sy,sz,now_crd)
    glog.append(nowdiff)
    glogdis.append(ldis)
    print str(n_iter)+'\t'+str(nowdiff)+'\t' +str(temp)
    print 'nowt:'+str(nowt)+'nowf:' + str(nowf)
    if temp == 0:
        break

job_server.destroy()

#(nowt,nowf,nowdiff,temp,step) = one_core(all_n_iter,nowt,nowf,nowdiff,temp,step,ltheta,lfai,ldiff,sx,sy,sz)

dur = time.time()-time1
print 'dur time: ' + str(dur)  + '\t' +str(n_iter)+'\t'+str(dur/n_iter)
print temp
print nowt
print nowf

#
plt.subplot(121)
plt.plot(range(len(glogdis)),glogdis)
plt.subplot(122)
plt.plot(range(len(glogdiff)),glogdiff)
