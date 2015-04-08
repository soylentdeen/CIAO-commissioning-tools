import scipy
import pyfits
import numpy
import matplotlib.pyplot as pyplot

datafile = '/diska/data/SPARTA/2015-04-07/TWHydra_55/TWHydra_55.fits'

data = pyfits.getdata(datafile)

slopes = data.field(4)
intensities = data.field(3)


flux = numpy.sum(intensities, axis=1)

xslopes = numpy.sum(slopes[:,0::2]*intensities, axis=1)
yslopes = numpy.sum(slopes[:,1::2]*intensities, axis=1)


fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

ax.plot(xslopes/flux, yslopes/flux)
ax.set_xlabel("Field Runout - Arcseconds")
ax.set_ylabel("Field Runout - Arcseconds")

fig.show()


