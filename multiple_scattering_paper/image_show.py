#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os
import re
import sys
import subprocess

# First wrote by Yun Zhao on 2013.10.24
# Last editted by Yun Zhao on 2013.11.4
# Yun Zhao is awesome!

##################### show image from a list #########################

def image_check():
    """ Show image from a list file with hdfsee """
    """ Last editted by Yun Zhao Oct 25"""
    filename = sys.argv[1]
    f = open(filename)
    for line in f:
        iamge_show = 'hdfsee '+line.strip()+' -g geom.geom -i 1000 -b 1'
        subprocess.call(iamge_show, shell=True)





##### read selected diffraction pattern with geometrical file #####

def read_image():
    """ Show the selected image with hdfsee """
    """ Input are the image list file and the order of image which you want to see"""
    """ Last editted by Yun Zhao on Oct 28 """
    filename =sys.argv[1]
    nth_image = int(sys.argv[2])
    f = open (filename)
    lines = f.readlines()
    f.close()
    print lines[nth_image]
    image_show='hdfsee '+lines[nth_image].strip()+' -g june24-v2.geom'
    print image_show
    subprocess.call(image_show,shell=True)



###############    Show predicted bragg spots #########################
def peak_prediction():
    """ Check the predicted spots in nth frame """
    """ Input file: stream file """
    """ Last editted by Yun Zhao at Oct 31,2013"""
#    filename = sys.argv[0]
    filename = 'good_data01.stream'
    f = open (filename)
    fo = open ('listfile.tmp','w')
    start = 0
    for line in f:
#        imagenamematch = re.search('Image filename: ',line.strip())
#        if imagenamematch:
#           imagepath = line.strip()
        if line.startswith('Image filename: '):
           imagepath = line[len('Image filename: '):]
           print imagepath
        if line.strip() == 'Peaks from peak search':
           start = 1
        if line.strip() == 'End of peak list':
           start = 0
           fo.close()
        if start == 1:
           fo.write(line)
        else:
           fo.close
        image_show='hdfsee '+imagepath.strip()+' -g june24-v2.geom'+' --peak-overlay=listfile.tmp'
        subprocess.call(image_show,shell=True)
        os.unlink('listfile.tmp')
        fo=open('listfile.tmp','w')
        start = 0



#####################     hkl file format    #########################
def full_hkl():
    """ convert the P1 hkl list to full reciprocal space """
    """ Last editted by Yun Zhao at Nov 6,2013 """
    filename = sys.argv[1]
    f1 = open(filename)
    f2 = open('full.hkl','w')
    for line in f1:
        f2.write(line)
        if re.search("[0-9]",line):      
           reflection=line.split()
           h =-int(reflection[0])
           k =-int(reflection[1])
           l =-int(reflection[2])
           intensity = float(reflection[3]) # Friedel pair share the same intensity
           phase =360-float(reflection[4]) # Friedel pair have the opposite phase
#          print h,"\n",k,"\n",l,"\n",intensity,"\n", phase
#          reflection = ' {0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(h,k,l,intensity,phase,0.0,1,0.0,0.0)
           reflection = "%(h)3i %(k)3i %(l)3i %(intensity)10.2f %(phase)8.2f %(sigma_I)10.2f %(counts)7i %(fs)6.1f %(ss)6.1f\n" % {"h":h,"k":k,"l":l,"intensity":intensity, "phase":phase, "sigma_I":0.0,"counts":1,"fs":0,"ss":0}
           f2.write(reflection)
    f1.close()
    f2.close()


def peak_check():
    """ Check those spots identified as a peak in imosflm """





def stat_plot():
    """ plot the statistics of indexing results (lattice constant and angles) """




    




