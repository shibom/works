#!/usr/local/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import math

#name ="/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0135-bg2/r0135-radialstack-detector0-class0-stack94.h5"

def resolution(rad):
    lamda = 1.30
    det = 108e-3
    angle = np.sin(0.5*np.arctan(rad/det))
    return (lamda/(2*angle))


def stack(filename):
    fd = h5py.File(filename)
    data = fd['data']
    data = data['data']
    data = np.array(data)
    radval = np.zeros(data.shape[1]);
    dist = []; count = [];
    for i in range(data.shape[1]):
        radval[i] += data[:,i].sum()
        radval[i] /= data.shape[0]
        dist.append(i*110e-6)
    return radval, dist

at_pix = []

infile = open('fullstack.txt', 'r')
for line in infile:
    line = line.split()
    radial, dis = stack(line[0])
    at_pix.append(radial[421])
    plt.plot(dis, radial)

infile.close()

plt.xlabel("Detector_radius(m)", fontsize=24, fontweight='bold')
plt.ylabel("radialAvg of Intensity", fontsize=24, fontweight='bold')

#plt.show()



plt.figure(2)
plt.plot(at_pix, '-o')
plt.xlabel("#stacks", fontsize=20, fontweight='bold')
plt.ylabel("Intensity_at_325pixel", fontsize=20, fontweight='bold')
plt.show()

at_pix = np.asarray(at_pix)
print np.std(at_pix, axis=0)

'''
for line in infile:
    line = line.split()
    radial, dis = stack(line[0])
    rad_at_400 = radial[300:500]
    fig.plot(rad_at_400)

plt.show()
'''
'''
for i in range(col.size):
    dist = (col[i]**2 + row[i]**2)**0.5
    dist *= 110e-6
    r.append(dist)   
#    rbin = round(dist)
#    sumval[rbin] += pos[i]
#    counter[rbin] += 1
#    if (counter[rbin] != 0):
#       sumval[rbin] /= counter[rbin]

plt.hist(r, 1000,  weights=pos, histtype='step')
plt.show()i
'''
