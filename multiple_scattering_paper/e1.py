#!/usr/bin/env python

import random

infile = open('junk.hkl', 'r')
ofile = open('error1k_tr.hkl', 'w')
ofile.write("  h   k   l          I    phase   sigma(F)  counts  fs/px  ss/px\n")
lst = []
all = infile.readlines()

for lines in all:
    lines = lines.split()
    h = int(lines[0])
    k = int(lines[1])
    l = int(lines[2])
    err = random.randint(-1000,1000)
    err = (100 + err)/100.0
    I = float(lines[3])
    I = I*err
   # F = I**0.5
    phase = float(lines[4])
    reflection = "%(h)3i %(k)3i %(l)3i %(I)10.2f %(phase)8.2f %(sigma_I)10.2f %(counts)7i %(fs)6.1f %(ss)6.1f\n" % {"h":h,"k":k,"l":l,"I":I, "phase":phase, "sigma_I":1.0,"counts":1,"fs":0,"ss":0}
    ofile.write(reflection)

ofile.write("End of reflections\n")

ofile.close()
infile.close()
