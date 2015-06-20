#!/usr/local/bin/env python

from psana import *
import matplotlib.pyplot as plt


file = open('TT_hist.txt', 'w')

ds = DataSource("exp=cxig0915:run=118")

dtype=Camera.FrameV1
src= Source('DetInfo(CxiDsu.0:Opal1000.0)')

time_delay = []; cal_list = [];
events = []; pix_pos = [];

for evtnum, evt in enumerate(ds.events()):
   #use "print evt.keys()" to see list of all types/sources in an event
   pos = ds.env().epicsStore().getPV('CXI:TTSPEC:FLTPOS').value(0)
   delay = ds.env().epicsStore().getPV('LAS:FS5:VIT:FS_TGT_TIME_DIAL').value(0)
   delay *=10e6;
   off = 5*pos + delay;
   pix_pos.append(delay); #cal_list.append(calib);
   time_delay.append(abs(off));
   file.write(str(abs(off)) + "\n") 
#   gdet = evt.get(dtype, src)
  # if gdet is None: continue
  # f.write(str(gdet.TimeTool.ConfigV1
  # delay = ds.env().epicsStore().getPV('LAS:FS5:VIT:FS_TGT_TIME_DIAL').value(0)
   #time_delay.append(delay) 
   events.append(evtnum)
   if evtnum > 2000:
      break
 # gdet = evt.
 # if gdet is None: continue
 # file.write(str(gdet.f_11_ENRC()) + "\n")
  #intensity_shot.append(gdet.f_11_ENRC())
 # events.append(evtnum)
 # if evtnum >= 2000:
 #    break
#f.close()
#print len(time_delay), len(pix_pos);

file.close()

plt.figure(1)
plt.plot(pix_pos)


plt.figure(2)
plt.hist(time_delay, 1000, histtype='step')
plt.xlabel('time_delay (fs)', fontsize=20, fontweight='bold')
plt.ylabel('Counts', fontsize=20, fontweight='bold')
#plt.savefig("signal_per_shot.png")
plt.show()

  #print "%6d:" % evtnum, id, gdet.f_11_ENRC()
 # if evtnum >= 10:
 #   break
