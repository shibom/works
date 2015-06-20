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

gasdata = np.loadtxt('gmd2val.txt')
meangas = gasdata.mean(0)
print meangas

name = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0135-bg2/r0135-radialstack-detector0-class0-stack12.h5"

def stack(filename):
    fd = h5py.File(filename)
    data = fd['data']
    data = data['data']
    data = np.array(data)
    datascaled = np.array(data)
    print data.shape
    sumval = np.zeros(data.shape[1])
    counter = np.zeros(sumval.size)
    dist = np.zeros(sumval.size)
    x = np.arange(0.,980.,1.)
    for j in range(gmd2sort.size-1):
        datascaled[j,:] /= gmd2sort[j]/meangas  # normalizing to the mean gmd value
    plt.plot(x[100:300],data[100:300,320],x[100:300],datascaled[100:300,320])
#    plt.plot(x,data[:,320],x,datascaled[:,320])
    plt.xlabel("#frames")
    plt.ylabel("Intensity_at_320pix")
    plt.show()

    print np.std(data[:,320])/np.mean(data[:,320])
    print np.std(datascaled[:,320])/np.mean(datascaled[:,320])

    for i in range(data.shape[0]):   # not shape[1] - you should be dividing by how many radials tehre are, not how many pixels.
       	#data[:,i] /= gmd2sort[i]
       	sumval[i] += data[:,i].sum()
       	dist[i] += i*110e-6
        sumval[i] /= data.shape[0]   # you were previously dividing by 1 repeatedly
    return sumval, dist

radial, dis = stack(name)
plt.figure(2)
plt.plot(radial)
plt.show()

