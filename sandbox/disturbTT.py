import scipy
import numpy
import matplotlib.pyplot as pyplot
import pyfits
import VLTTools

ciao = VLTTools.VLTConnection(simulate=False)
ciao.disturbTT()
#print "This is a test"
