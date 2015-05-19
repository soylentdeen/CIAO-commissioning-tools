import scipy
import matplotlib.pyplot as pyplot
import numpy
import pyfits
import os

datadir = '/diska/data/SPARTA/2015-04-17/TTMap_1/'

pupilShiftFile = open(os.path.expanduser('~')+'/data/TTMapping/TTMCommandSequence.txt', 'r')

Tip = []
Tilt = []

gradX = []
gradY = []
"""
for line in pupilShiftFile:
    l = line.split()
    Tip.append(float(l[0]))
    Tilt.append(float(l[1]))
    LoopName = 'TTMapLoop_'+str(Tip[-1])+'_'+str(Tilt[-1])+'/TTMapLoop_'+str(Tip[-1])+'_'+str(Tilt[-1])+'.fits'
    loopData = pyfits.getdata(datadir+LoopName)
    PixelName = 'TTMapPixels_'+str(Tip[-1])+'_'+str(Tilt[-1])+'/TTMapPixels_'+str(Tip[-1])+'_'+str(Tilt[-1])+'.fits'
    PixelData = pyfits.getdata(datadir+PixelName)
    images = pixelData.field(3)

    gradients = loopData.field(4)
    gradX.append(numpy.average(gradients[:,0::2], axis=0))
    gradY.append(numpy.average(gradients[:,1::2], axis=0))

gradX = numpy.array(gradX)
gradY = numpy.array(gradY)
Tip = numpy.array(Tip)
Tilt = numpy.array(Tilt)
"""

PixelName = 'TTMapPixels_'+str(-0.02)+'_'+str(0.0)+'/TTMapPixels_'+str(-0.02)+'_'+str(0.0)+'.fits'
pixelData = pyfits.getdata(datadir+PixelName)
images = pixelData.field(3)
avg = numpy.average(images, axis=0)
diff = images[0]-images[-1]
indexes = range(100)
std = numpy.zeros(100)
for i in indexes:
    std[i] = numpy.std(images[i] - images[-1])


fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#ax.matshow(avg.reshape(72,72))
#ax.matshow(diff.reshape(72,72))
ax.plot(indexes, std)
#order = numpy.argsort(Tip)
#ax.plot(Tip[order], gradX[:,0][order])
#ax.scatter(gradX[:,0], gradY[:,0])
#ax.scatter(Tip, Tilt)

#for i in range(len(gradX[0])):
#    ax.plot(Tip[order], gradX[:,i][order])


fig.show()
