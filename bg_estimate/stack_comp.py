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
    return ((2*angle)/lamda)


def stack(filename):
    fd = h5py.File(filename)
    data = fd['data']
    data = data['data']
    data = np.array(data)
    radval = np.zeros(data.shape[1]);
    sigval = np.zeros(data.shape[1]);
    dist = []; count = [];
    for j in range(data.shape[0]):
        check = data[j,:].sum()/data.shape[1]
        if (check > 10):   # check and exclude blanks with mean intensity of 30;
           count.append(j)
    for i in range(data.shape[1]):
        for j in range(len(count)):
            radval[i] += data[count[j],i]
        radval[i] /= len(count)
        sigval[i] += np.var(data[0:len(count),i]) # calculate variance for each column in each stack
        dist.append(i*110e-6)
    dist = np.array(dist)
    reso = resolution(dist)
    return radval, sigval, reso

flist = ['fullstack.txt', 'fullstack16_aga.txt']
rad_list = []
sig_list = []

fig = plt.figure(1)
ax1 = fig.add_subplot(111)

for file in flist:
    infile = open(file, 'r')
    all_rads = np.zeros(1181)
    all_vars = np.zeros(1181)
    count = 0
    for line in infile:
        line = line.split()
        radial, sigval, dis = stack(line[0])
        all_rads[:] += radial
        all_vars[:] += sigval
        count += 1
    infile.close()
    all_rads /= count   # take avg. for radials from each stack
    all_vars /= count   # take avg. over the variance from each stack
    all_sigs = np.sqrt(all_vars) # calculate std. dev. from variance
    rad_list.append(all_rads)
    sig_list.append(all_sigs)

all_rad_lcp = rad_list[0]
all_rad_aga = rad_list[1]
all_sig_lcp = sig_list[0]
all_sig_aga = sig_list[1]    

scale = all_rad_aga[750]/all_rad_lcp[750] # scale factor for  the agarose radials at 2-A resolution

all_rad_aga /= scale
all_sig_aga /= (scale)**0.5

print np.max(all_sig_aga)
print np.max(all_rad_aga)

#sigma_lcp = np.std(all_rad_lcp)
#sigma_aga = np.std(all_rad_aga)

ax1.plot(dis, all_rad_lcp)
ax1.plot(dis, all_rad_aga)
ax1.fill_between(dis, all_rad_lcp - all_sig_lcp, all_rad_lcp + all_sig_lcp, color='gray', alpha=0.5)
ax1.fill_between(dis, all_rad_aga - all_sig_aga, all_rad_aga + all_sig_aga, color='gray', alpha=0.5)

ax1.set_xlabel("1/d(=2sin$\Theta$/$\lambda$) (1/$\AA$)", fontsize=20, fontweight='bold')
ax1.set_ylabel("radialAvg of Intensity", fontsize=20, fontweight='bold')
ax1.legend(['LCP', 'Agarose'])
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


plt.show()


'''
plt.figure(2)
plt.plot(at_pix, '-o')
plt.xlabel("#stacks", fontsize=20, fontweight='bold')
plt.ylabel("Intensity_at_325pixel", fontsize=20, fontweight='bold')
plt.show()
'''


