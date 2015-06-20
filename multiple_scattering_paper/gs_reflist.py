#!/usr/bin/env python

import random
import csv

infile = open('gs_reflist2.csv', 'rb')
ofile = open('gs2.hkl', 'w')
ofile.write("  h   k   l          I    phase   sigma(F)  counts  fs/px  ss/px\n")
lst = []
all = csv.reader(infile)

for lines in all:
    #lines = lines.split()
    h = int(lines[1])
    k = int(lines[2])
    l = int(lines[3])
   # I = float(lines[3])
    ra = random.randint(0, 3000)
    I = ra
    phase = '-'
    reflection = "%(h)3i %(k)3i %(l)3i %(I)10.2f %(phase)8.2s %(sigma_I)10.2f %(counts)7i %(fs)6.1f %(ss)6.1f\n" % {"h":h,"k":k,"l":l,"I":I, "phase":'-', "sigma_I":100.0,"counts":1,"fs":0,"ss":0}
    ofile.write(reflection)

ofile.write("End of reflections\n")

ofile.close()
infile.close()
   

