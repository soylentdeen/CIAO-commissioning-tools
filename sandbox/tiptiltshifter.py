import numpy as np
import pyfits
import sys
import matplotlib.pyplot as pyplot

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

fnum = sys.argv[1]

path = '/diska/data/SPARTA/30March2015/BetaPic_'+fnum+'/'
filename='BetaPic_'+fnum+'.fits'
a = pyfits.getdata( path + filename )
fl = a.field(3)
vTT = (a.field(6))[:,0]

subF = np.array([1,6,7,13,14,15,22,23,24,31,32,33,39,40,41,48,49,50,57,58,64]) - 1
sub0 = np.array([2,8,16,25,34,42,51,59,65]) - 1


flF = np.average(  fl[:, subF], axis=1 )
fl0 = np.average(  fl[:, sub0], axis=1 )

ttmin = np.min(vTT)

flF_min = np.average(flF[ np.where( np.abs(vTT-ttmin)<0.001 )[0] ])
flF_max = np.average(flF[ np.where( np.abs(vTT-ttmin)>0.001 )[0] ])

fl0_min = np.average(fl0[ np.where( np.abs(vTT-ttmin)<0.001 )[0] ])
fl0_max = np.average(fl0[ np.where( np.abs(vTT-ttmin)>0.001 )[0] ])

print flF_max - flF_min
print fl0_max - fl0_min

ax.plot(flF)
ax.plot(fl0)
fig.show()

