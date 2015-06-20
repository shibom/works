#!/usr/local/bin env python

import h5py
import numpy as np
import matplotlib.pyplot as plt

file = h5py.File('../geometry/g0915_mpi.h5', 'r')
x = file['x']
y = file['y']
z = file['z']
x = np.array(x); y = np.array(y); z = np.array(z);

rx = np.square(x)
ry = np.square(y)
rz = np.square(z)

r = np.sqrt(rx + ry + rz)/110e-6

mask = np.ones(r.shape)

mask[r>900] = 0
plt.imshow(mask)
plt.show()

f = h5py.File('shroud_900pix.h5', 'w')
data = f.create_group("data")
data.create_dataset("data", data=mask)
f.close()

