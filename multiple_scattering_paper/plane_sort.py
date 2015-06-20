#!/usr/bin/env python

ifile = open('p1_F.hkl', 'r')
ofile = open('sort_p1.hkl', 'w')

all = ifile.readlines()
ofile.write("  h   k   l          F    phase   sigma(F)  counts  fs/px  ss/px\n")
for lines in all:
    lines = lines.split()
    h = int(lines[0])
    k = int(lines[1])
    l = int(lines[2])
    if (l == 0):
       for n in range(-10,10):
           if (h == 2*n):
              for t in range(-10,10):
                  if (k == 2*t):
                      F = float(lines[3])
                      #F = I**0.5
                      phase = float(lines[4])
                      refl_even = "%(h)3i %(k)3i %(l)3i %(F)10.2f %(phase)8.2f %(sigma_I)10.2f %(counts)7i %(fs)6.1f %(ss)6.1f\n" % {"h":h,"k":k,"l":l,"F":F, "phase":phase, "sigma_I":0.0,"counts":1,"fs":0,"ss":0}
                      ofile.write(refl_even)
           elif (h == (2*n+1)):
              for t in range(-10,10):
                  if (k == (2*t+1)):
                     F = float(lines[3])
                     #F = I**0.5
                     phase = float(lines[4])
                    # print h, k, l, F, phase
                     refl_odd = "%(h)3i %(k)3i %(l)3i %(F)10.2f %(phase)8.2f %(sigma_I)10.2f %(counts)7i %(fs)6.1f %(ss)6.1f\n" % {"h":h,"k":k,"l":l,"F":F, "phase":phase, "sigma_I":0.0,"counts":1,"fs":0,"ss":0}
                     ofile.write(refl_odd)
 
