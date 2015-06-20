#!/usr/bin/env python

import csv

infile = open('gs_reflist2.csv', 'rb')

all = csv.reader(infile)
#print all
h = []; k = []; l = []
for lines in all:
    print lines[3]
    #lines = lines.split()
    #print lines[0]
#    h.append(int(lines[2]))
#    k.append(int(lines[3]))
#    l.append(int(lines[4]))
     
#print h,k,l
    #h = lines[2]
    #k = lines[3]
    #l = lines[4]
    #print lst
