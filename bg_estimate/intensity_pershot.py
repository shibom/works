#!/usr/local/bin/env python

from psana import *
import matplotlib.pyplot as plt


file = open('pershote79_198k.txt', 'w')

ds = DataSource("exp=cxie7914:run=135")

dtype=Bld.BldDataFEEGasDetEnergyV1
src= Source('BldInfo(FEEGasDetEnergy)')

intensity_shot = [];
events = []
for evtnum, evt in enumerate(ds.events()):
  # use "print evt.keys()" to see list of all types/sources in an event
  id = evt.get(EventId)
  gdet = evt.get(dtype,src)
  if gdet is None: continue
  file.write(str(gdet.f_11_ENRC()) + "\n")
  #intensity_shot.append(gdet.f_11_ENRC())
  events.append(evtnum)
  if evtnum >= 198000:
     break
f.close()
'''
plt.plot(events, intensity_shot)
plt.savefig("signal_per_shot.png")
plt.show()
'''
  #print "%6d:" % evtnum, id, gdet.f_11_ENRC()
 # if evtnum >= 10:
 #   break
