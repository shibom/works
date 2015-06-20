#!/usr/bin/env python

import os, sys
import matplotlib.pyplot as plt

file = open(sys.argv[1], 'r')
inline = []; rwork = []; rfree = []; step = []; res = [] 

if (sys.argv[2] == 'step'):
   for line in file:
       if (line.strip() == "REMARK  stage     r-work r-free  xray_target_w  xray_target_t"):
           break
   for line in file:
        if (line.strip() == "REMARK ------------------------------------------------------------------------"):
            break
        inline.append(line)
   for line in inline:
       line = line.split()
       rwork.append(line[2])
       rfree.append(line[3])
   rwork = rwork[1:]
   rfree = rfree[1:]
   plt.plot(rwork, '-o', linewidth=2)
   plt.plot(rfree, 'r-o', linewidth=2)
   plt.xlabel('refinement steps', fontsize=24, fontweight='bold')
   plt.legend(['$R_{work}$', '$R_{free}$'])
   plt.show()

elif (sys.argv[2] == 'reso'):
     for line in file:
         if (line.strip() == "REMARK   3   BIN  RESOLUTION RANGE  COMPL.    NWORK NFREE   RWORK  RFREE"):
             break
     for line in file:
         if (line.strip() == "REMARK   3"):   
            break
         inline.append(line)
     for line in inline:
         line = line.split()
         rwork.append(float(line[9]))
         rfree.append(float(line[10]))
         res.append(float(line[5]))

     plt.plot(res, rwork, '-o', linewidth=2)
     plt.plot(res, rfree, 'r-o', linewidth=2)
     plt.gca().invert_xaxis()
     plt.xlabel('Resolution(Ang)', fontsize=24, fontweight='bold')
     plt.ylabel('R-factors', fontsize=24, fontweight='bold')
     plt.legend(['$R_{work}$', '$R_{free}$'], loc='upper left')
     plt.show()
else:
     exit()
file.close()

