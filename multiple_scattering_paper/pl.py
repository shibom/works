#!/usr/bin/env python

ifile = open('all_refl.hkl', 'r')
ofile = open('sort.hkl', 'w')

all = ifile.readlines()

for lines in all:
    lines = lines.split()
    h = int(lines[0])
    k = int(lines[1])
    l = int(lines[2])
    for n in range(-10, 10):
        if (h == 2*n):
           for t in range(-10, 10):
               k = 2*t
               F
