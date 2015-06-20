#!/usr/local/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt

name ="/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0135-bg2/r0135-radialstack-detector0-class0-stack8.h5"
#namie ="/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0135-bg2/r0135-detector0-class0-sum.h5"

fd = h5py.File(name)
data = fd['data']
data = data['data']
data = np.array(data)
#pos = data[data > 0]

#loc = np.where(data > 0)
#row = loc[0]; col = loc[1];
sumval = np.zeros(data.shape[1])
#sumval = []
r = 0; val = []; c = [];
for j in range(data.shape[0]):
    check = data[j,:].sum()/data.shape[1]
    if (check > 20):
        c.append(j); 

print len(c)
for i in range(data.shape[1]):         
    for j in range(len(c)):
        sumval[i] += data[c[j],i]
    sumval[i] /= len(c)
#print r
#for i in range(len(val)):
#    sumval.append(val[i]/r) 
'''
for i in range(data.shape[1]):
    sumval.append(data[:,i].sum()/r)
          # sumval[i] /= data.shape[0]
'''
print len(sumval)
'''
val = []
for i in range(len(r)):
    val.append(sumval[i]/len(r))
'''
'''
for i in range(data.shape[1]):
    sumval[i] += data[:,i].sum()
    sumval[i] /= data.shape[0]
'''
plt.plot(sumval)
plt.show()

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
