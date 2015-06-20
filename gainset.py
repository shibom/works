#!/usr/local/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt

file = h5py.File('LG09_Schmidt_150318-fudged.h5', 'r');
data = file['data'];
data = data['data'];
tim = np.array(data);
print np.max(tim);

tim[tim > 6.8] = 6.82;

plt.imshow(tim);
plt.show()

f = h5py.File('correctgain_LG09.h5', 'w');
data = f.create_group("data");
data.create_dataset("data",data=tim);
f.close()

