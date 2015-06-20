#!/usr/bin/env python


import sys

infile = open(sys.argv[1], 'r')
ofile = open(sys.argv[2], 'w')

lst = []
all = infile.readlines()

for line in all:
     line = line.split()
     lst.append(line)

for ii in range(1, 101):
     #print lst[20*ii]
      ll = lst[20*ii] 
      ll = ",".join(ll).replace(",","\t")
      s = str(ll);
      ofile.write(s + "\n")

infile.close()
ofile.close()
