#!/usr/local/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import math

#name ="/reg/d/psdm/cxi/cxie7914/scratch/shibom/cheetah/hdf5/r0135-bg2/r0135-radialstack-detector0-class0-stack94.h5"

name = "/reg/d/psdm/cxi/cxig4914/scratch/shibom/cheetah/hdf5/r0016-bg/r0016-radialstack-detector0-class0-stack8.h5";

pix_loc = 421;

def resolution(rad):
    lamda = 1.30
    det = 108e-3
    angle = np.sin(0.5*np.arctan(rad/det))
    return (lamda/(2*angle))


def diffuse_intensity(filename):
    fd = h5py.File(filename)
    data = fd['data']
    data = data['data']
    data = np.array(data)
    radval = [];
    dist = []; count = [];
    for j in range(data.shape[0]):
        check = data[j,:].sum()/data.shape[1]
        if (check > 30):
           count.append(j)
           diffuse = data[j, pix_loc]
           radval.append(diffuse)
    return radval,count


infile = open('fullstack16_aga.txt', 'r')
#no = 0;
plt.figure(1)
all_radials = []; all_count = [];

for line in infile:
    line = line.split()
 #   no += 1;
  #  plt.figure(no)
    radial, count = diffuse_intensity(line[0])
    all_radials.append(radial)
    all_count.append(count)
    plt.scatter(count, radial)
    plt.xlim(0,len(count))

infile.close()

#radial, count = diffuse_intensity(name)
#plt.scatter(count, radial)

plt.xlabel("#frames", fontsize=24, fontweight='bold')
plt.ylabel("peak_Intensity_at_421pix", fontsize=24, fontweight='bold')
#plt.xlim(0,len(count))

all_radials_inone_list = [item for sublist in all_radials for item in sublist]

all_radials_inone_list = np.array(all_radials_inone_list)

histo_at_diffuse, binEdges = np.histogram(all_radials_inone_list, bins=100)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])

sigma = np.std(histo_at_diffuse);
print sigma;

plt.figure(2)

plt.plot(bincenters, histo_at_diffuse, '-')
plt.fill_between(bincenters, histo_at_diffuse - sigma, histo_at_diffuse + sigma, color='green')
plt.xlabel("peak_intensity_at_421_pixel", fontsize=24, fontweight='bold')
plt.ylabel("Frequency", fontsize=24, fontweight='bold')


plt.show()


'''
plt.figure(2)
plt.plot(at_pix, '-o')
plt.xlabel("#stacks", fontsize=20, fontweight='bold')
plt.ylabel("Intensity_at_325pixel", fontsize=20, fontweight='bold')
plt.show()
'''
#at_pix = np.asarray(at_pix)

#print np.std(at_pix, axis=0)

'''
for line in infile:
    line = line.split()
    radial, dis = stack(line[0])
    rad_at_400 = radial[300:500]
    fig.plot(rad_at_400)

plt.show()
'''
'''
for i in range(col.size):
    dist = (col[i]**2 + row[i]**2)**0.5
    dist *= 110e-6
    r.append(dist)   
#    rbin = round(dist)
#    sumval[rbin] += pos[i]
#    counter[rbin] += 1
#    if (counter[rbin] != 0):
#       sumval[rbin] /= counter[rbin]

plt.hist(r, 1000,  weights=pos, histtype='step')
plt.show()i
'''
