#!/usr/local/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py

#name_aga = "/reg/d/psdm/cxi/cxig4914/scratch/shibom/cheetah/hdf5/r0016-sb2/r0016-detector0-class0-sum.h5"
#name_lcp = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0135-bg2/r0135-detector0-class0-sum.h5"
name_lcp_mask = "lcp135_masked.h5"

file = h5py.File(name_lcp_mask)
data = file['data']
#data = data['radialAverage']
data = data['data']
data = np.array(data)

p_val = data[data > 0]
loc = np.where(data > 0)
row = loc[1]; col = loc[0];
#print loc
#print p_val
#[row, col] = np.indices(data.shape)
#print row.shape, col.shape
cntX = 740; cntY = 760; 
pix_size = 110e-6;
dist = np.sqrt(np.square(row-cntX) + np.square(col-cntY))
row = loc[1]; col = loc[0];
radial = np.zeros(p_val.shape)
counter = np.zeros(p_val.shape)
#r = np.zeros(p_val.shape)
print row
print col
#radial = radial.flatten()
#counter = counter.flatten()
#data = data.flatten()

for i in range(dist.size):
    rbin = round(dist[i])
    radial[rbin] += p_val[i]
    counter[rbin] += 1

'''
for i in range(row.size):
    for j in range(col.size):
        distX = ((col[j] - cntX))**2
        distY = ((row[i] - cntY))**2
        dist = (distX + distY)**0.5
    r[i] += dist
print r
'''
'''   
    rbin = round(dist)
    radial[rbin] += p_val[i]
    counter[rbin] += 1
#    r.append(rbin)
'''
for i in range(radial.size):
    if (counter[i] != 0):
       radial[i] /= counter[i]

#print data.shape
#print r


#r = np.array(r)
#print r.shape

#hist, binedge = np.histogram(r, bins=1000, weights=radial)
#print hist
#plt.step(binedge[:-1], hist)
#plt.ylim(10e+2, 10e+8)'''

plt.figure(2)
plt.hist(dist, bins=1000, weights=radial, histtype='step')


plt.show()

