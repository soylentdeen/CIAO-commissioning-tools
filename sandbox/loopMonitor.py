import scipy
import pyfits
import numpy
import matplotlib.pyplot as pylot

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
s2z = pyfits.getdata('S2Z.fits')

