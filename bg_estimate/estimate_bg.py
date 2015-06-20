#!/usr/local/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import h5py

def resolution(rad):
    lamda = 1.30
    det = 108e-3
    angle = np.sin(0.5*np.arctan(rad/det))
    return (lamda/(2*angle))


#name_lcp = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0183-bg2/r0183-radialstack-detector0-class0-stack6.h5"
#name_aga = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0189-bg1/r0189-radialstack-detector0-class0-stack6.h5"

name_aga = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0137-bg1/r0137-detector0-class0-sum.h5"
name_lcp = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0183-bg2/r0183-detector0-class0-sum.h5"
name_lcp_mask = "lcp135_masked.h5"
#name_aga = "/reg/d/psdm/cxi/cxig4914/scratch/shibom/cheetah/hdf5/r0016-sb2/r0016-detector0-class0-sum.h5"

cntX_aga = 740; cntY_aga = 776; pix_size = 110e-6;
cntX_lcp = 886; cntY_lcp = 886;

def estimator(filename, cnt1, cnt2):
    file = h5py.File(filename, 'r')
    data = file['data']
    data = data['data']
    data = np.array(data)
    p_val = data[data > 0]
    loc = np.where(data > 0)
    row = loc[0]; col = loc[1];
    r = []
    for i in range(col.size):
        distX = ((col[i] - cnt1)*pix_size)**2
        distY = ((row[i] - cnt2)*pix_size)**2
        dist = (distX + distY)**0.5
        r.append(dist)
    return p_val, r 



aga_val, r_aga = estimator(name_aga, cntX_aga, cntY_aga); 
lcp_val, r_lcp = estimator(name_lcp_mask, cntX_lcp, cntY_lcp);

r_aga = np.array(r_aga)
#reso_aga = resolution(r_aga)

r_lcp = np.array(r_lcp)
#reso_lcp = resolution(r_lcp)
'''
plt.figure(1)
plt.hist(reso_aga, bins=1000,  normed=True, weights=aga_val/len(aga_val), histtype='step')
plt.hist(reso_lcp, bins=1000,  normed=True, weights=lcp_val/len(lcp_val), histtype='step')
plt.gca().invert_xaxis()
plt.xlabel("Resolution (A)", fontsize=24, fontweight='bold')
plt.ylabel("Avg. Intensity[counts]", fontsize=24, fontweight='bold')
plt.legend(['Agarose', 'LCP'], loc='upper left')
#bc = 0.5*(nbins[1:]+nbins[:-1])
#plt.plot(bc,events, '-')
plt.show()
'''
fig = plt.figure(1)
ax1 = fig.add_subplot(111)
n1, bin1, patch1 = ax1.hist(r_aga, bins=1000, weights=aga_val, histtype='step')
n2, bin2, patch2 = ax1.hist(r_lcp, bins=1000, weights=lcp_val, histtype='step')
ax1.set_xlabel("detector_radius (m)", fontsize=24, fontweight='bold')
ax1.set_ylabel("Avg. Intensity[counts]", fontsize=24, fontweight='bold')
plt.legend(['Agarose', 'LCP'], loc='upper left')

ax2 = ax1.twiny()
ax2.set_xlabel("Resolution (A)", fontsize=24, fontweight='bold')
ax2.set_xticklabels(['Inf','7.11','3.68','2.60','2.07','1.78','1.59','1.47'])

plt.show()
