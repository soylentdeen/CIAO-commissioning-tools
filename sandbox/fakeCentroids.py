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
            if LA.SLapertureMap[i][j]:
                stamp = frame[xstart:xstart+8, ystart:ystart+8]
                xcent = numpy.sum(
                      numpy.sum(stamp, axis=0)*xweight)/numpy.sum(stamp)
                ycent = numpy.sum(
                      numpy.sum(stamp, axis=1)*yweight)/numpy.sum(stamp)
                retval.append(xcent)#+xstart+4)
                retval.append(ycent)#+ystart+4)
    return retval

    

#Pixels = pyfits.getdata('squarePSF.fits')
Pixels = pyfits.getdata('TipTiltMotion.fits')

centroids = []

for pix in Pixels:
    centroids.append(extractCentroids(pix))
#centroids.append(extractCentroids(Pixels[0]))

centroids = numpy.array(centroids)


#ax.scatter(mn[0::2], st[0::2], color='b', marker = '+', label='X')
#ax.scatter(mn[1::2], st[1::2], color='g', marker = 'o', label='Y')
for i in range(len(centroids[0])/2):
    ax.plot(centroids[1:,i*2]-centroids[0,i*2], centroids[1:,i*2+1]-centroids[0,i*2+1], color= 'b', marker = '+', label='Mean Error')
#ax.scatter(st[0::2], st[1::2], color = 'g', marker = 'o', label='Standard Deviation')

ax.set_title("Measured Centroids")
ax.set_xlabel("X (Pixels)")
ax.set_ylabel("Y (Pixels)")
#ax.set_xbound(-5e-8, 1e-7)
#ax.set_ybound(-5e-8, 1e-7)
#ax.legend()
fig.show()
fig.savefig("ExtractedCentroids.png")
