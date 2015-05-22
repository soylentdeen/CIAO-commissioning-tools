import scipy
import matplotlib.pyplot as pyplot
import numpy
import pyfits
from NGCTools import lensletArray

LA = lensletArray()

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

def extractCentroids(frame):
    retval = []
    xweight = numpy.arange(8)-3.5
    yweight = numpy.arange(8)-3.5
    for i in range(9):
        xstart = i*8
        for j in range(9):
            ystart = j*8
            if LA.apertureMap[i][j]:
                stamp = frame[xstart:xstart+8, ystart:ystart+8]
                xcent = numpy.sum(
                      numpy.sum(stamp, axis=0)*xweight)/numpy.sum(stamp)
                ycent = numpy.sum(
                      numpy.sum(stamp, axis=1)*yweight)/numpy.sum(stamp)
                retval.append(xcent)
                retval.append(ycent)
    return retval

    

datadir = '/diska/data/SPARTA/2015-05-22/'

LoopData = pyfits.getdata(datadir+'LoopData_23/LoopData_23.fits')
PixelData = pyfits.getdata(datadir+'PixelData_5/PixelData_5.fits')


Centroids = LoopData.field(4)
LoopFrames = LoopData.field(0)

PixelFrames = PixelData.field(0)
Pixels = PixelData.field(-1)

commonPix = numpy.in1d(PixelFrames, LoopFrames)
commonLoop = numpy.in1d(LoopFrames, PixelFrames)

LoopFrames = LoopFrames[commonLoop]
Centroids = Centroids[commonLoop]
print len(LoopFrames)

PixelFrames = PixelFrames[commonPix]
Pixels = Pixels[commonPix]
print len(PixelFrames)

errors = []

for cent, pix in zip(Centroids, Pixels):
    errors.append(extractCentroids(pix.reshape(72,72))-cent)

errors = numpy.array(errors)

mn = numpy.mean(errors, axis=0)
st = numpy.std(errors, axis=0)

#ax.scatter(mn[0::2], st[0::2], color='b', marker = '+', label='X')
#ax.scatter(mn[1::2], st[1::2], color='g', marker = 'o', label='Y')
ax.scatter(mn[0::2], mn[1::2], color= 'b', marker = '+', label='Mean Error')
ax.scatter(st[0::2], st[1::2], color = 'g', marker = 'o', label='Standard Deviation')

ax.set_title("Linear Centroiding Estimation Algorithm Testing")
ax.set_xlabel("X difference CoG - LCE (Pixels)")
ax.set_ylabel("Y difference CoG - LCE (Pixels)")
#ax.set_xbound(-5e-8, 1e-7)
#ax.set_ybound(-5e-8, 1e-7)
ax.legend()
fig.show()
fig.savefig("LCE_Test.png")
