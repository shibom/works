#!/usr/bin/env python

ifile1 = open("s1.hkl", 'r')
ifile2 = open('s3.hkl', 'r')

s1 = {}; s3 = {}
h1 = []; h2 = []; I1 = []; I2 = []
Is1 = []; Is2 = []
all1 = ifile1.readlines()
all2 = ifile2.readlines()

for lines in all1:
    lines = lines.split()
    h1.append(lines[0:3])
    I1.append(lines[3:])

for lines in all2:
    lines = lines.split()
    h2.append(lines[0:3])
    I2.append(lines[3:])
#print h1, h2
h2 = h2[0:1000]

for i in range(len(h2)):
    for j in range(len(h2)):
        if (h1[i] == h2[j]):        
           Is1.append(h1[i]+I1[i])
           Is2.append(h2[j]+I2[j])
           del h1[i]
           del h2[j]
           print len(h2) 
ofile1 = open('js1.hkl', 'w')
ofile2 = open('js3.hkl', 'w')

ofile1.write("  h   k   l          I    phase   sigma(I)  counts  fs/px  ss/px\n")
ofile2.write("  h   k   l          I    phase   sigma(I)  counts  fs/px  ss/px\n")
for el in Is1:
    h = int(el[0])
    k = int(el[1])
    l = int(el[2])
    I = float(el[3])
    sig = float(el[5])
    reflection = "%(h)3i %(k)3i %(l)3i %(I)10.2f %(phase)8.2s %(sigma_I)10.2f %(counts)7i %(fs)6.1f %(ss)6.1f\n" % {"h":h    ,"k":k,"l":l,"I":I, "phase":'-', "sigma_I":sig,"counts":1,"fs":0,"ss":0}
    ofile1.write(reflection)

ofile1.close()

for el in Is2:
    h = int(el[0])
    k = int(el[1])
    l = int(el[2])
    I = float(el[3])
    sig = float(el[5])
    reflection = "%(h)3i %(k)3i %(l)3i %(I)10.2f %(phase)8.2s %(sigma_I)10.2f %(counts)7i %(fs)6.1f %(ss)6.1f\n" % {"h":h    ,"k":k,"l":l,"I":I, "phase":'-', "sigma_I":sig,"counts":1,"fs":0,"ss":0}
    ofile2.write(reflection)

ofile2.close()
#for lines in all1:
#    (key, val) = lines.split()
#    s1[int(key[0:3])] = val[3:]
#print s1


