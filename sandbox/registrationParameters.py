import pyfits
import matplotlib.pyplot as pyplot
import sys
import numpy
import termios
import tty
import time


#CMFile = "/diska/home/ciaomgr/data/Recn.REC1.CM.fits"
#CMFile = "/diska/home/ciaomgr/data/controlMatrix.fits"
#CM = pyfits.getdata(CMFile)
#TT2HOfile = "/diska/home/ciaomgr/data/HOCtr.TT_TO_HO.fits"
#TT2HO = pyfits.getdata(TT2HOfile)
#loopFile = "/diska/data/SPARTA/30March2015/"+sys.argv[1]+'/'+sys.argv[1]+'.fits'
loopFile = "/home/deen/Data/GRAVITY/13March2015/Beaumont_9/Beaumont_9.fits"
loopData = pyfits.getdata(loopFile)

t = loopData.field(1)
gradients = loopData.field(4)
hodm = loopData.field(5)
ttm = loopData.field(6)
actuators = numpy.append(hodm, ttm, axis=1)
st = numpy.std(hodm, axis=0)
fig = pyplot.figure(1)
fig.clear()
fig.show()
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.4])
ax2 = fig.add_axes([0.1, 0.5, 0.8, 0.4])
ActPos1, = ax1.plot(actuators[0,:], c='r')
Gradx, = ax2.plot(gradients[0,0::2])
Grady, = ax2.plot(gradients[0,1::2])


#ax1 = fig.add_axes([0.1, 0.1, 0.4, 0.8])
#ax2 = fig.add_axes([0.5, 0.1, 0.4, 0.8])

#ax1.matshow(hodm, aspect='auto')
#ax2.plot(st)
#ax2.plot(hodm[:,5])

#fig.show()
