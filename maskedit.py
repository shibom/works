#!/usr/bin/env python

import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as pimg
mask = np.ones((1480,1552))
mask = np.uint16(mask)




#h5FileName="r0184-detector0-class1-sum.h5"
h5FileName="shroud_mask.h5"

f = h5py.File(h5FileName, "r")
data = f["data"]
tim = data["data"]
tim = np.array(tim)
sh = tim.copy()
f.close()

sh[925:1014,584:775] = 1
sh[740:924,388:640] = 1
sh[370:739,776:880] = 1 #q2a4 and q2a6
sh[555:739,388:470] = 0 #q1a6
sh[1110:1175,388:581] = 1 #q1a12
'''
#dark[np.where(dark < minthresh)] = minthresh;
#dark[np.where(dark > maxthresh)] = maxthresh;
#mask[1110:1294,0:193] = 0  # q0a12 and q0a13
#mask[1110:1294,194:387] = 0 # q0a13
mask[370:554,388:470] = 0  # q1a4 and q1a5
mask[555:739,388:470] = 0  # q1a6 and q1a7
#mask[740:924,388:581] = 0  # q1a8
#mask[830:924,582:650] = 0 # 1/4th of q1a9
mask[925:1109,388:583] = 0  # q1a10 and q1a11
mask[1014:1109,584:775] = 0 # q1a11
#mask[1110:1174,388:581] = 0 # q1a12
#mask[370:739,776:830] = 0 # 1/3rd of q2a4 and q2a6
mask[925:1109,776:973] = 0 # q2a10 and q2a11
mask[1049:1109,974:1163] = 0 # 1/3rd of q2a11
mask[740:923,0:95] = 0 # 1/4th of q0a8
mask[925:1109, 0:387] = 0  # full of q0a10 and q0a11
mask[925:1109,873:960] = 0
#mask[970:1109,0:48] = 0 # 1/4th of q0a10
#mask[1090:1109,49:193] = 0
mask[925:1109,1164:1357] = 0 # q3a10
mask[1110:1205,1164:1357] = 0 #half of q3a12
mask[740:924,776:800] = 0
mask[800:924,801:860] = 0
#mask[880:924,861:969] = 0


imgplot = plt.imshow(mask, cmap=pimg.cm.RdBu, interpolation="nearest",norm=None)
#imgplot.set_clim(100.0,25000.0)
plt.colorbar()
plt.show()
'''
print("writing file...\n")
f = h5py.File("shroud-edit.h5","w")
data = f.create_group("data");
data.create_dataset("data",data=sh);
f.close()

