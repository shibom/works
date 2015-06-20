#!/usr/local/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import h5py
#import math

def resolution(rad):
    lamda = 1.33
    det = 118e-3
    angle = np.sin(0.5*np.arctan(rad/det))
    return (lamda/(2*angle))

#name1 = ["/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0189-bg1/r0189-radialstack-detector0-class0-stack0.h5", "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0189-bg1/r0189-radialstack-detector0-class0-stack1.h5", "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0189-bg1/r0189-radialstack-detector0-class0-stack2.h5", "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0189-bg1/r0189-radialstack-detector0-class0-stack3.h5", "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0189-bg1/r0189-radialstack-detector0-class0-stack4.h5", "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0189-bg1/r0189-radialstack-detector0-class0-stack5.h5", "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0189-bg1/r0189-radialstack-detector0-class0-stack6.h5", "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0189-bg1/r0189-radialstack-detector0-class0-stack7.h5", "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0189-bg1/r0189-radialstack-detector0-class0-stack8.h5", "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0189-bg1/r0189-radialstack-detector0-class0-stack9.h5", "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0189-bg1/r0189-radialstack-detector0-class0-stack10.h5"]

#name_lcp = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0183-bg2/r0183-radialstack-detector0-class0-stack6.h5"
#name_aga = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0189-bg1/r0189-radialstack-detector0-class0-stack6.h5"

name_aga = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0189-bg1/r0189-detector0-class0-sum.h5"
name_lcp = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0183-bg2/r0183-detector0-class0-sum.h5"

def readfile(filename):
    file = h5py.File(filename, 'r')
    data = file['data']
    data = data['data']
    data = np.array(data)
    p_val = data[data > 0]
    loc = np.where(data > 0)
    return p_val, loc

#p_val = p_val.tolist()
aga_val, aga_loc = readfile(name_aga)
lcp_val, lcp_loc = readfile(name_lcp)

#loc = np.where(data > 0)
aga_row = aga_loc[0]; lcp_row = lcp_loc[0];
aga_col = aga_loc[1]; lcp_col = lcp_loc[1];

r_lcp = []; r_aga = [];
#plt.plot(data.sum(axis=1))
#plt.show()

cntX = 885; cntY = 885; pix_size = 110e-6;

for i in range(aga_row.size):
    distX = ((aga_col[i] - cntX)*pix_size)**2
    distY = ((aga_row[i] - cntY)*pix_size)**2
    dist = (distX + distY)**0.5
    r_aga.append(dist)

for i in range(lcp_row.size):
    distX = ((lcp_col[i] - cntX)*pix_size)**2
    distY = ((lcp_row[i] - cntY)*pix_size)**2
    dist = (distX + distY)**0.5
    r_lcp.append(dist)

r_aga = np.array(r_aga)
reso_aga = resolution(r_aga)
r_lcp = np.array(r_lcp)
reso_lcp = resolution(r_lcp)

plt.hist(reso_aga, bins=1000, range=(1.8,6), normed=True, weights=aga_val/len(aga_val), histtype='step')
plt.hist(reso_lcp, bins=1000, range=(1.8,6), normed=True, weights=lcp_val/len(lcp_val), histtype='step')
plt.gca().invert_xaxis()
plt.xlabel("Resolution (A)", fontsize=24, fontweight='bold')
plt.ylabel("Avg. Intensity[counts]", fontsize=24, fontweight='bold')
plt.legend(['Agarose', 'LCP'], loc='upper left')
#bc = 0.5*(nbins[1:]+nbins[:-1])
#plt.plot(bc,events, '-')

plt.show()
