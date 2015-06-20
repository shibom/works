#!/usr/local/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import h5py

name ="/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0135-bg2/r0135-radialstack-detector0-class0-stack76.h5"
name1 = "radials/edit1_radstack74.h5"

file = h5py.File(name)
data = file['data']
data = data['data']
tim = np.array(data)
file.close()
 
plt.figure(1)
plt.imshow(tim, interpolation='nearest', norm=None)
plt.xlabel("Radius (pixels)", fontsize=24, fontweight='bold')
plt.ylabel("#frames",fontsize=24, fontweight='bold')
plt.colorbar()
plt.show()

print tim.shape
#junk = np.zeros(tim.shape)
row = [];
for i in range(tim.shape[0]):
    check = tim[i,:].sum()/tim.shape[1]
    if (check < 20):
       row.append(i)

print row[222], row[221]
count = row[0]
for ii in range(1, len(row)+1):
#    print tim.shape
    tim = np.delete(tim, count, 0)
    count = row[ii] - 1


#for i in range(600,970):
#    tim = np.delete(tim, i, 0)

#for i in range(0,140):
#    tim = np.delete(tim, i, 0)

print tim.shape

plt.figure(2)
plt.imshow(tim, interpolation='nearest', norm=None)
plt.xlabel("Radius (pixels)", fontsize=24, fontweight='bold')
plt.ylabel("#frames",fontsize=24, fontweight='bold')
plt.colorbar()
plt.show() 

f = h5py.File("radials/test76.h5","w")
data = f.create_group("data");
data.create_dataset("data",data=tim);
f.close()

