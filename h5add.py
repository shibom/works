#!/usr/bin/env python

import argparse
import h5py
import os, sys
import numpy as np
import matplotlib.pyplot as plt

def readh5(filename):
    if os.stat(filename).st_size > 0:
       f = h5py.File(filename, "r")
       if f["data"] is None:
          print "junks"
          pass
       data = f["data"]
       tim = data["rawdata"]
       tim = np.array(tim)
       peaks = f['processing']
       peaks = peaks['cheetah']
       peaks = peaks['peakinfo-assembled']
       peaks = np.array(peaks)
       im = tim + peaks
    else:
       print "couldn't open file \n"
    f.close()
    return im

def writeh5(filename,mydat):
	f = h5py.File(filename,"w")
	data = f.create_group("data");
	data.create_dataset("rawdata",data=mydat);
	f.close()
	return
'''
parser = argparse.ArgumentParser()
parser.add_argument("-i", action="append", dest="h5FileNames", type=str, nargs='+', help="intput files")
parser.add_argument("-o", action="store", dest="h5SavePath", type=str, nargs=1, help="output file")
args = parser.parse_args()

h5FileNames = args.h5FileNames[0]
h5FileName = h5FileNames[0]
h5SavePath = args.h5SavePath[0]

mask = readh5(h5FileName)
'''
h5list = [];
powder = np.zeros((1480,1552))

infile = open('h5list_r16.txt', 'r')
for lines in infile:
    lines = lines.split()
    h5list.append(lines[0])

   
for i in range(1,2000):
	data = readh5(h5list[i])
	powder += data
	
writeh5('mypowder_tr1.h5', powder)

