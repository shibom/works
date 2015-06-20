#!/usr/local/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import math


#name = "/reg/d/psdm/cxi/cxig4914/scratch/shibom/cheetah/hdf5/r0016-bg/r0016-radialstack-detector0-class0-stack8.h5";
file_list = ['fullstack.txt', 'fullstack16_aga.txt'];
pix_loc_lcp = 323;
pix_loc_aga = 421;

def resolution(rad):
    lamda = 1.30
    det = 108e-3
    angle = np.sin(0.5*np.arctan(rad/det))
    return (lamda/(2*angle))


def diffuse_intensity(filename, pix_loc):
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

for file in file_list:
    infile = open(file, 'r')
    all_radials = []; all_count = [];
    for line in infile:
        line = line.split()
        if (file == 'fullstack.txt'):
           radial, count = diffuse_intensity(line[0],pix_loc_lcp)
        else: 
           radial, count = diffuse_intensity(line[0],pix_loc_aga)
        all_radials.append(radial)
        all_count.append(count)

    all_radials_inone_list = [item for sublist in all_radials for item in sublist]

    all_radials_inone_list = np.array(all_radials_inone_list)

    histo_at_diffuse, binEdges = np.histogram(all_radials_inone_list, bins=100)
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    sigma = np.std(histo_at_diffuse)
    plt.plot(bincenters, histo_at_diffuse, '-')
    plt.fill_between(bincenters,histo_at_diffuse - sigma, histo_at_diffuse+sigma)
    infile.close()


plt.xlabel("peak_intensity_at_diffuse_ring", fontsize=20, fontweight='bold')
plt.ylabel("Counts", fontsize=20, fontweight='bold')
plt.legend(['LCPring_at_4.5-A', 'Agarose_ring_at_3.0-A'], loc='upper right')

plt.show()

