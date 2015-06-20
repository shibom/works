#!/usr/local/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py

name_aga = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0135-bg2/r0135-detector0-class0-sum.h5"

file = h5py.File(name_aga, 'r')
data = file['data']
data = data['data']
tim = np.array(data)

mask = np.ones(data.shape)

mask[300:850,1450:1500] = 0
mask[1040:1060,1700:1720] = 0
mask[1290:1320,360:380] = 0

plt.imshow(tim*mask, interpolation='nearest', norm=None)
plt.show()

f = h5py.File("lcp135_masked.h5","w")
data = f.create_group("data");
data.create_dataset("data",data=mask*tim);
f.close()
