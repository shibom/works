#!/usr/local/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

data = np.loadtxt("pershote79_198k.txt")
plt.plot(data)
plt.xlabel("#Events", fontsize=24, fontweight='bold')
plt.ylabel("Intensity/shot", fontsize=24, fontweight='bold')
plt.show()

