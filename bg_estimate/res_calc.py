#!/usr/local/bin/env python

import math
import sys

def resolution(rad):
    lamda = 1.30
    det = 108e-3
    angle = math.sin(0.5*math.atan(rad/det))
    return (lamda/(2*angle))

dist = [0.02,0.04,0.046,0.05,0.06,0.08,0.10,0.12,0.14]

for r in dist:
    print resolution(r)

