#!/usr/bin/env python

#This handy Cell histogram script is written by Sbasu

import os, sys
import numpy as np
from pylab import *
import matplotlib.pyplot as plt

def stdev(vals):
    avg = mean(vals)
    diff = [(val - avg)**2 for val in vals]
    return mean(diff)**0.5

def mode(thelist):
  counts = {}
  for item in thelist:
    counts [item] = counts.get (item, 0) + 1
  maxcount = 0
  maxitem = None
  for k, v in counts.items ():
    if v > maxcount:
      maxitem = k
      maxcount = v
  return maxitem

infile = open(sys.argv[1], 'r')
all_lines = infile.readlines()
a = []
b = []
c = []
al = []
be = []
ga = []

for lines in all_lines:
    lines = lines.split()
    a.append(float(lines[0]))
    b.append(float(lines[1]))
    c.append(float(lines[2]))
    al.append(float(lines[3]))
    be.append(float(lines[4]))
    ga.append(float(lines[5]))



fig = plt.figure()
ax = fig.add_subplot(2,3,1)
ax.hist(a, 100, histtype='stepfilled')
ax.set_title('mode=%f, sigma=%f' %(mode(a), stdev(a)))

ax = fig.add_subplot(2,3,2)
ax.hist(b, 100, histtype='stepfilled')
ax.set_title('mode=%f, sigma=%f' %(mode(b), stdev(b)))

ax = fig.add_subplot(2,3,3)
ax.hist(c, 100, histtype='stepfilled')
ax.set_title('mode=%f, sigma=%f' %(mode(c), stdev(c)))

ax = fig.add_subplot(2,3,4)
ax.hist(al, 100, histtype='stepfilled')
ax.set_title('mode=%f, sigma=%f' %(mode(al), stdev(al)))

ax = fig.add_subplot(2,3,5)
ax.hist(be, 100, histtype='stepfilled')
ax.set_title('mode=%f, sigma=%f' %(mode(be), stdev(be)))

ax = fig.add_subplot(2,3,6)
ax.hist(ga, 100, histtype='stepfilled')
ax.set_title('mode=%f, sigma=%f' %(mode(ga), stdev(ga)))
     
plt.show()

infile.close()

