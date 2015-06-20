#!/usr/bin/env python

import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

def readh5(filename):
   f = h5py.File(filename, "r")
   data = f["data"]
   tim = data["data"]
   tim = np.array(tim)
   im = tim.copy()
   f.close()
   return im

def writeh5(filename,mydat):
   f = h5py.File(filename,"w")
   data = f.create_group("data");
   data.create_dataset("data",data=mydat);
   f.close()
   return
mask = readh5('cxig0915-r0010-badpixels-combined.h5')
mask *= readh5('shroud-edit.h5')
#mask *= readh5('ringmask.h5')


#for i in range(1,2):
 #   im = readh5(sys.argv[i])
  #  mask *= im

writeh5("cxig0915-r0010-badpixels-combined-shroud.h5",mask)

