#!/usr/bin/env python

#my handy script for plotting everything! - Sb

import sys, os
from pylab import *
import pylab as p
import numpy as np


infile = open(sys.argv[1], 'r')
all_lines = infile.readlines()
metric = []
reso = []
for lines in all_lines:
    lines = lines.split()
    metric.append(lines[1])
    reso.append(lines[3])

metric = np.array(metric[1:])
reso = np.array(reso[1:])
metric = metric.astype(np.float32)
#metric = metric/100.0
reso = reso.astype(np.float32)

p.plot(reso, metric, '-o', linewidth=2)
p.gca().invert_xaxis()
p.xlabel('Resolution(Ang)', fontsize=24, fontweight='bold')
p.ylabel('Rsplits(%)', fontsize=24, fontweight='bold')
p.show()
infile.close()
