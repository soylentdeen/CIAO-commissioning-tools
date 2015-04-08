import scipy
import numpy
import matplotlib.pyplot as pyplot
import pyfits
import VLTTools

ciao = VLTTools.VLTConnection(simulate=False)
ciao.disturbHO(disturbType="SINE", rng=0.5, max=0.95, period=20.0, actNum=50)
#print "This is a test"
