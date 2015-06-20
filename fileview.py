#!/usr/bin/env python

import sys,os
import h5py
import numpy as np
import matplotlib.pyplot as plt
#from PSCalib.GeometryAccess import img_from_pixel_arrays

#fname = "lg09.h5"

def readh5(fname, dataset):
    fh = h5py.File(fname, 'r')
    data = fh[dataset]
    data = np.array(data)

    return data

def load_geom(pixmap):
    fh = h5py.File(pixmap, 'r')
    pix_size = 110e-6
    x = fh['/x']
    y = fh['/y']
    z = fh['/z']
    x = np.array(x); y = np.array(y);
    xmin = x.min() - pix_size/2; 
    ymin = y.min() - pix_size/2;
    iX = np.array((x - xmin)/pix_size, dtype=np.uint); 
    iY = np.array((y - ymin)/pix_size, dtype=np.uint);

    return iX, iY

def map_to_pixmap(data,iX,iY):
   
    iX = iX.flatten()
    iY = iY.flatten() 
    
    xs = iX.max() + 1
    ys = iY.max() + 1
    verbose = 1
    data = data.flatten()

    im = verbose*np.ones((xs,ys), dtype=np.float32)
    im[iX, iY] = data

    return im

def write_with_geom(outname, image):
   
    f = h5py.File(outname,"w")
    data = f.create_group("data");
    data.create_dataset("data",data=image);

     
fname = sys.argv[1]
pixmap = sys.argv[2]

data = readh5(fname, '/data/data')
[iX, iY] = load_geom(pixmap)
image = map_to_pixmap(data, iX, iY)

print "writing paneled h5 files\n"
outname = sys.argv[3]

write_with_geom(outname, image)

plt.imshow(image)
plt.show()

