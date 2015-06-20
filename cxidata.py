#!/usr/local/bin/env python

from psana import *
import numpy as np
import matplotlib.pyplot as plt


#cfg = "cspad_psana_module.cfg"
cfg="cxig0915.cfg"

setConfigFile(cfg)

ds = DataSource('exp=cxig0915:run=130:idx')
run = ds.runs().next()
time = run.times()

dtype = CsPad.DataV2
src = Source('DetInfo(CxiDs1.0:Cspad.0)')


for tm in time:
    evt = run.event(tm)
    '''
    cspad = evt.get(dtype, src)
    a = []
    for i in range(0,4):
        quad = cspad.quads(i)
        d = quad.data()
        a.append(np.vstack([d[i] for i in range(0,8)]))
 
    frame_raw = np.hstack(a)
    '''
    frame_reconstructed = evt.get(ndarray_float64_2,src,'reconstructed')
    
    if frame_reconstructed is None: continue
    print frame_reconstructed.shape

    frame_calib =  evt.get(ndarray_float64_3,src,'calibrated')
    
    if frame_calib is None: continue

    frame_calib = frame_calib.reshape((4,8*185,388))
    frame_stack = np.hstack([frame_calib[i,:,:] for i in range(0,4)])

    plt.imshow(frame_stack)

#    plt.imshow(frame_reconstructed)
    plt.show()
    

    


