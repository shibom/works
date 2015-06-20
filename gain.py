#!/usr/local/bin env python

import h5py
import numpy as np
import matplotlib.pyplot as plt

file = h5py.File('../geometry/cspad-yef.h5', 'r')
x = file['x']
y = file['y']
z = file['z']
x = np.array(x); y = np.array(y); z = np.array(z);

rx = np.square(x)
ry = np.square(y)
rz = np.square(z)

r = np.sqrt(rx + ry + rz)/110e-6

gain = np.ones(r.shape)

gain[r<240] = 6.3
plt.imshow(gain)
plt.show()

f = h5py.File('gain_sb.h5', 'w')
data = f.create_group("data")
data.create_dataset("data", data=gain)
f.close()

