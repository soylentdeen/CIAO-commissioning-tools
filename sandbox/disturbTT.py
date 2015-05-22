import scipy
import numpy
import matplotlib.pyplot as pyplot
import pyfits
import VLTTools
import sys

datadir = '/diska/data/SPARTA/2015-05-22/'

ciao = VLTTools.VLTConnection(simulate=False, datapath=datadir)
ciao.disturbTT(tip=float(sys.argv[1]), tilt=float(sys.argv[2]), waveshape='SINE')
#print "This is a test"
