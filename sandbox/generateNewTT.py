import scipy
import pyfits
import numpy
import VLTTools
import sys

ciao = VLTTools.VLTConnection(simulate=False, datapath = '/diska/data/SPARTA/2015-04-08/')

#ciao.averageTTPositions(sys.argv[1])
ciao.measureNewTTRefPositions("Junk")
