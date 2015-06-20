#!/usr/local/bin/env python

import numpy as np
import matplotlib.pyplot as plt 
import h5py

name_aga = "/reg/d/psdm/cxi/cxig4914/scratch/shibom/cheetah/hdf5/r0016-sb2/r0016-detector0-class0-sum.h5"
name_lcp = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0135-bg2/r0135-detector0-class0-sum.h5"
#name_lcp_mask = "lcp135_masked.h5"

def readfile(name):
    file = h5py.File(name)
    data = file['data']
    data = data['radialAverage']
    #data = data['data']
    data = np.array(data)
    p_val = data[data > 0]
    loc = np.where(data > 0)
    rad = loc[0]; dist = np.zeros(rad.size);
    for i in range(rad.size):
        dist[i] += rad[i]*110e-6
    return dist, p_val

dist_aga, aga_val = readfile(name_aga)
dist_lcp, lcp_val = readfile(name_lcp)
#x = [0,0.02,0.04,0.06,0.08,0.10,0.12,0.14]
#y = [0,2e5,4,
x = range(0,8)
y = range(0,7)

fig = plt.figure(1)
ax1 = fig.add_subplot(111)
ax1.plot(dist_aga, aga_val)
ax1.plot(dist_lcp, lcp_val)
#n1, bin1, patch1 = ax1.hist(r_aga, bins=1000, weights=aga_val, histtype='step')
#n2, bin2, patch2 = ax1.hist(r_lcp, bins=1000, weights=lcp_val, histtype='step')

ax1.set_xlabel("detector_radius (m)", fontsize=24, fontweight='bold')
ax1.set_ylabel("Avg. Intensity[counts]", fontsize=24, fontweight='bold')

locs,labels = plt.xticks()
plt.xticks(locs, map(lambda x: "%3.2f" % x, locs))
plt.xlim(0,0.14)
locs,labels = plt.yticks()
plt.yticks(locs, map(lambda x: "%de5" % x, locs*1e-5))
plt.tick_params(axis="y", labelsize=15)
plt.tick_params(axis="x", labelsize=15)
plt.legend(['Agarose', 'LCP'], loc='upper right')

ax2 = ax1.twiny()
ax2.set_xlabel("Resolution (A)", fontsize=24, fontweight='bold')
ax2.set_xticklabels(['Inf','7.11','3.68','2.60','2.07','1.78','1.59','1.47'])
#plt.plot(dist, p_val)
plt.show()
