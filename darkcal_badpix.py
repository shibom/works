#!/usr/bin/env python

import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt

mask = np.ones((1480,1552))
mask = np.uint16(mask)



# Let's mask out the regions of the dark image that have very low and very high dark current.  Low dark current probably means a dead pixel.  High dark current probably means the pixel is hot.

h5FileName="badpixelmap.h5"
f = h5py.File(h5FileName, "r")
data = f["data"]
tim = data["data"]
tim = np.array(tim)
dark = tim.copy()
f.close()

mask[771:780,310:318] = 0
mask[754:764,902:908] = 0
mask[929:937,708:712] = 0
mask[1140:1144,318:323] = 0
#mask[dark > 5000] = 0

#mask[sigma > 50] = 0

#plt.ion()
plt.imshow(mask*dark,interpolation="nearest",norm=None)
plt.show()

print("writing file...\n")
f = h5py.File("badpix_141018.h5","w")
data = f.create_group("data");
data.create_dataset("data",data=mask*dark);
f.close()

