#!/usr/local/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py


frame = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0135-bg2/frames.txt"

fr_info = np.loadtxt(frame, skiprows=1, delimiter=',', usecols=(7,8))
gmd1 = fr_info[:,0]; gmd2 = fr_info[:,1];
gmd2sort = gmd2[12000:12980]
f = open("gmd2val.txt", 'w')
for i in range(gmd2sort.shape[0]):
    f.write(str(gmd2sort[i]) + '\n')

f.close()

name = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0135-bg2/r0135-radialstack-detector0-class0-stack12.h5"

def stack(filename):
    fd = h5py.File(filename)
    data = fd['data']
    data = data['data']
    data = np.array(data)
    print data.shape
    sumval = np.zeros(data.shape[1])
    counter = np.zeros(sumval.size)
    dist = np.zeros(sumval.size)
   # for j in range(gmd2sort.size):
    #    data[j,:] /= gmd2sort[j]
    plt.plot(data[:,320])
    plt.xlabel("#frames")
    plt.ylabel("Intensity_at_320pix")
    plt.show()
    for i in range(data.shape[1]):
        #data[:,i] /= gmd2sort[i]
        sumval[i] += data[:,i].sum()
        counter[i] += 1
        dist[i] += i*110e-6
        if (counter[i] != 0): 
           sumval[i] /= counter[i]
    return sumval, dist

radial, dis = stack(name)
plt.figure(2)
plt.plot(radial)
plt.show()

