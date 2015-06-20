#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

infile = open('stat.dat', 'r')

all = infile.readlines()
Rw = []; Rf = []; F = []; delRw = []; delRf = []
lst = [];
for lines in all:
    lines = lines.split()
    F.append(lines[1]); Rw.append(lines[3]); Rf.append(lines[6])
    delRw.append(lines[4]); delRf.append(lines[7])
print Rf
F = F[1:]; Rw = Rw[1:]; Rf = Rf[1:]; delRw = delRw[1:]; delRf = delRf[1:]

plt.figure(0)
plt.plot(F, Rw, '-o', linewidth=2)
plt.plot(F, Rf, 'r-o', linewidth=2)
plt.xlabel('Error in |F|(%)', fontsize=24, fontweight='bold')
plt.ylabel('R-factors', fontsize=24, fontweight='bold')
plt.legend(['$R_{work}$', '$R_{free}$'], loc='upper left')
plt.show()

plt.figure(1)

plt.plot(F, delRw, '-o', linewidth=2)
plt.plot(F, delRf, 'r-o', linewidth=2)
plt.xlabel('Error in |F|(%)', fontsize=24, fontweight='bold')
plt.ylabel('$\delta R-factors$(%)', fontsize=24, fontweight='bold')
plt.legend(['$\delta R_{work}$', '$\delta R_{free}$'], loc='upper left')
plt.show()

