#!/usr/local/bin/env python

import sys, os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import math

# Let's define some useful constants first..

det = 118e-3           # 118-mm detector distance
lamda = 1.3e-10        #1.3-A X-ray wavelength
pix_size = 110e-6      #110-um each pixel size

# function to calculate Bragg angle..

def cal_angle(reso):
    oneoverd = (1.0/reso)*1e10 # resolution using m^-1 unit
    bragg = (lamda*oneoverd)/2 # lamda/2d = sin(theta) Bragg's Law
    return math.asin(bragg)

#function to calculate radii in pixel units of resolution rings..

def cal_radius(angle):
    angle = 2*angle                          # scattering angle is 2 times of Bragg angle!
    npix = (det*math.tan(angle))/pix_size      
    npix = round(npix)
    return npix

#function to draw circles..

def circle(r,phi,h,k):
    return h+r*np.cos(phi), k+r*np.sin(phi)

#function to make annulus masks..

def annulus_mask(rmin, rmax,h,k):
    #ymin,xmin = np.ogrid[-h+rmin:h+rmin+1, -k+rmin:k+rmin+1]
    #ymax,xmax = np.ogrid[-h+rmax:h+rmax+1, -k+rmax:k+rmax+1]
    #anal = (xmax**2 - xmin**2) + (ymax**2 - ymin**2) <= (rmax**2 - rmin**2)
    y1,x1 = np.ogrid[-h:1772-h, -k:1772-k]
    y2,x2 = np.ogrid[-h:1772-h, -k:1772-k]
    anal1 = (x1**2 + y1**2) < rmin**2    
    anal2 = (x2**2 + y2**2) < rmax**2
    anal = anal2 - anal1
    return anal

# Read the h5 file first..

#filename1 = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0137-bg1/r0137-detector0-class0-sum.h5"

filename1 = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0189-bg1/r0189-detector0-class1-sum.h5"
#filename2 = "/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0183-bg2/r0183-detector0-class0-sum.h5"

def read_data(file):
    fil = h5py.File(file, 'r')
    data = fil['data']
    data = data['data']
    tim = np.array(data)
    return tim

aga_bg = read_data(filename1)
aga_bg /= 10681

#lcp_bg = read_data(filename2)
#lcp_bg /= 12036
 
# lists for radii and angles of the circles 
r = []; theta = [];
reso = [10,8,5,4,3,2]    # i want to bin everything in those resolution bins..

for res in reso:
    scat = cal_angle(res)
    theta.append(2*scat)
    rad = cal_radius(scat)
    r.append(rad)
print r

# Let's define the origin of the image/circles

h = 885; k = 885;
#h = 739; k = 739;

BG_aga = []     # a list for background intensity
BG_lcp = []

kernel = np.zeros(aga_bg.shape)
#y,x = np.ogrid[-h:1772-h, -k:1772-k]
#mask = x**2 + y**2 <= 140**2
#mask1 = x**2 + y**2 <= 170**2
#m = mask1 - mask
#kernel[m] = 1
mask = annulus_mask(r[0],r[1],h,k)
kernel[mask] = 1
print np.mean(aga_bg*kernel)
aga_bg[mask] = 0
#avg_I = np.mean(tim*kernel)
#BG.append(avg_I)
'''
for i in range(-1,5):
    kernel = np.zeros(aga_bg.shape)
    mask1 = annulus_mask(r[i-1],r[i+1])
    kernel[mask1] = 1
    BG_aga.append(np.mean(aga_bg*kernel))
    BG_lcp.append(np.mean(lcp_bg*kernel))

plt.figure(1)
plt.plot(reso, BG_aga, 'o-',linewidth=2)
#plt.plot(reso, BG_lcp, 'o-r',linewidth=2)
plt.gca().invert_xaxis()
plt.xlabel("Resolution", fontsize=20, fontweight='bold')
plt.ylabel("mean background intensity", fontsize=20, fontweight='bold')
#plt.legend(['Agarose', 'LCP'], loc='upper right')
plt.show()
'''

plt.figure(2)
plt.imshow(aga_bg, interpolation='nearest', norm=None)
for rad in r:
    angle = np.arange(0,6.28,0.01)
    plt.plot(*circle(rad,angle,h,k), color='r')
plt.xlim(0,1772)
plt.ylim(1772,0)
plt.show()

