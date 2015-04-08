import pyfits
import matplotlib.pyplot as pyplot
import sys
import numpy

datafile1= "/diska/data/SPARTA/30March2015/"+sys.argv[1]+'/'+sys.argv[1]+'.fits'
datafile2= "/diska/data/SPARTA/30March2015/"+sys.argv[2]+'/'+sys.argv[2]+'.fits'

data1 = pyfits.getdata(datafile1)
data2 = pyfits.getdata(datafile2)

hodm1 = data1.field(5)
hodm2 = data2.field(5)
#st = numpy.std(hodm, axis=0)
print "Max Value - No SatMan: ", numpy.max(numpy.abs(hodm1))
print "Max Value - With SatMan: ", numpy.max(numpy.abs(hodm2))
fig = pyplot.figure(0)
fig.clear()
ax1 = fig.add_axes([0.1, 0.1, 0.4, 0.4])
ax2 = fig.add_axes([0.5, 0.1, 0.4, 0.4])
ax3 = fig.add_axes([0.1, 0.5, 0.4, 0.4])
ax4 = fig.add_axes([0.5, 0.5, 0.4, 0.4])

ax1.matshow(hodm1)
ax2.matshow(hodm2)
#ax2.plot(st)
ax3.plot(hodm1[:,50])
ax3.plot(hodm2[:,50])

fig.show()
