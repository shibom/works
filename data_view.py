import sys, os
import psdata
import numpy as np
import matplotlib.pyplot as plt
import time
import argparse


def calc_data(quads, quad_size):
    data = np.empty((4,8,185,388))

    for ii in xrange(quad_size[0]):
        q = quads(ii).data()
    for jj in xrange(8):
        data[ii,jj] = q[jj]
    return data

def plotting(data):
    plt.ion()
    plt.show(False)

    fig = plt.figure()

    for i in xrange(data.shape[0]):
        ax = fig.add_subplot(2,2,i)
        ax.imshow(data[i,4,:,:])
        ax.set_title("quadrant-%i" %i)
        plt.draw()

parser = argparse.ArgumentParser()
parser.add_argument("-e", "--exp", type=str, 
                        help='Experiment number')
parser.add_argument("-r", "--run", type=int, default=-1, 
                        help='Run number')
parser.add_argument("-nevt", "--nevents", type=int, default=2,
                        help='event number')

args = parser.parse_args()

cxi = psdata.psdata(exp=args.exp, run=args.run)

for k in xrange(args.nevents):
    cxi.next_event()
    quads = cxi.DscCsPad.evtData.CsPadData.quads
    quad_size = cxi.DscCsPad.evtData.CsPadData.quads_shape
    data_array = calc_data(quads, quad_size)
    plotting(data_array)
#    plt.pause(0.05)
    time.sleep(0.5)

'''for ii in xrange(quad_size[0]):
    q = quads(ii).data()
    for jj in xrange(8):
        data[ii,jj] = q[jj]

fig = plt.figure()

for i in xrange(data.shape[0]):
    ax = fig.add_subplot(2,2,i)
    ax.imshow(data[i,4,:,:])
    ax.set_title("quadrant-%i" %i)

plt.show()
'''
