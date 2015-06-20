#!/usr/bin/env python

ifile = open('slc.hkl', 'r')
ofile = open('tr2.hkl', 'w')

all = ifile.readlines()

for lines in all:
    lines = lines.split()
    h = int(lines[0])
    k = int(lines[1])
    l = int(lines[2])
    F = float(lines[3])
    h = k; k1 = -int(lines[0]); l = -l 
    phase = float(lines[4])
    refl1 = "%(h)3i %(k)3i %(l)3i %(F)10.2f %(phase)8.2f %(sigma_I)10.2f %(counts)7i%(fs)6.1f %(ss)6.1f\n" % {"h":h,"k":k1,"l":l,"F":F, "phase":phase, "sigma_I":0.0,"counts":1,"fs":0,"ss":0}
    ofile.write(refl1)
 
