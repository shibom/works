#!/usr/local/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import math

#name ="/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0135-bg2/r0135-radialstack-detector0-class0-stack94.h5"

def resolution(rad):
    lamda = 1.30
    det = 301.16e-3
    angle = np.sin(0.5*np.arctan(rad/det))
    return ((2*angle)/lamda)


def stack(filename):
    fd = h5py.File(filename)
    data = fd['data']
    data = data['data']
    data = np.array(data)
    radval = np.zeros(data.shape[1]);
    dist = []; count = []; 
    frm_wt_blank = data.shape[0]
    for j in range(data.shape[0]):
        check = data[j,:].sum()/data.shape[1]
        if (check > 20):
           count.append(j)
    for i in range(data.shape[1]):
        for j in range(len(count)):
            radval[i] += data[count[j],i]
        radval[i] /= len(count)
        dist.append(i*110e-6)
    dist = np.array(dist)
    reso = resolution(dist)
    return radval, reso, frm_wt_blank, len(count)

at_pix = []
all_rads = np.zeros(1173)
count = 0; nstack = 0; nframe = 0;

fig = plt.figure(1)
ax1 = fig.add_subplot(111)

infile = open('r377radial.txt', 'r')
for line in infile:
    line = line.split()
    radial, dis, stacks, frms = stack(line[0])
    all_rads[:] += radial
    count += 1
    nstack += stacks
    nframe += frms
    at_pix.append(radial[323])
    ax1.plot(dis, radial)

infile.close()
print nstack, nframe

ax1.set_xlabel("1/d(=2sin$\Theta$/$\lambda$) (1/$\AA$)", fontsize=20, fontweight='bold')
ax1.set_ylabel("radialAvg of Intensity", fontsize=20, fontweight='bold')
ax2 = ax1.twiny()

ax2xs = []
for x in ax1.get_xticks():
    tmp = 1/x
    num = "%3.2f" % tmp
    ax2xs.append(num)
ax2.set_xticks(ax1.get_xticks())
ax2.set_xbound(ax1.get_xbound())
ax2.set_xticklabels(ax2xs)
ax2.set_xlabel("Resolution ($\AA$)", fontsize=20, fontweight='bold')

all_rads /= count
sigma = np.std(all_rads)

plt.figure(2)
all_rads_one_array = np.vstack(all_rads)
histo, binEdges = np.histogram(all_rads_one_array, bins=100) 
bincenters = 0.5*(binEdges[1:] + binEdges[:-1])

plt.plot(bincenters, histo, '-')


fig1 = plt.figure(3)
ax3 = fig1.add_subplot(111)
ax3.plot(dis, all_rads)
ax3.fill_between(dis, all_rads - sigma, all_rads + sigma, color='gray')
ax3xs = ax3.get_xticks()

ax3.set_xlabel("1/d(=2sin$\Theta$/$\lambda$) (1/$\AA$)", fontsize=20, fontweight='bold')
ax3.set_ylabel("radialAvg of Intensity", fontsize=20, fontweight='bold')
ax4 = ax3.twiny()
ax4.set_xlabel("Resolution ($\AA$)", fontsize=20, fontweight='bold')
ax4.set_xticks(ax3xs)
ax4.set_xbound(ax3.get_xbound())
ax4.set_xticklabels(ax2xs)


plt.show()

'''
plt.figure(2)
plt.plot(at_pix, '-o')
plt.xlabel("#stacks", fontsize=20, fontweight='bold')
plt.ylabel("Intensity_at_325pixel", fontsize=20, fontweight='bold')
plt.show()
'''


